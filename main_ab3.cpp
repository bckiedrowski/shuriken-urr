#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <vector>

#include "ace.hpp"
#include "cell.hpp"
#include "constants.hpp"
#include "dosEstimator.hpp"
#include "dosUREstimator.hpp"
#include "estimator.hpp"
#include "estimator.hpp"
#include "fixedSourceCalc.hpp"
#include "kEigenvalueCalc.hpp"
#include "material.hpp"
#include "neutronData.hpp"
#include "particle.hpp"
#include "point.hpp"
#include "random.hpp"
#include "surface.hpp"
#include "urr.hpp"

#include "expIntegral.hpp"

int main() {

  util::RNG::ptr rng = std::make_shared<util::Rand>();
  const auto num_histories = 1e6;
  const auto num_batches   = 10;

  //cross sections and materials
  std::vector< std::shared_ptr<xs::aceData> > ace;

  ace.push_back( std::make_shared< xs::aceData >( "data/40090.710nc" ) );

  auto niso = 0;
  std::map< uint32_t, std::shared_ptr<xs::neutronData> > xs;
  {
  for ( auto a : ace ) {
    xs.emplace( a->zaid(), std::make_shared< xs::neutronData >( *a, niso ));
    ++niso;
  }
  ace.clear();
  }
  xs::xsBuffer xs_buffer( niso );

  //create slab geometry
  mc::material zr90;
  zr90.add( 4.3675e-2, xs[40090] );

  const auto plane_left  = std::make_shared< geom::plane >( "plane_left",  0, 1.0, 0.0, 0.0,  0.0 );
  const auto plane_right = std::make_shared< geom::plane >( "plane_right", 1, 1.0, 0.0, 0.0, 10.0 );

  const auto surfaces = std::vector< std::shared_ptr<geom::surface> >( { plane_left, plane_right } );

  const auto slab = std::make_shared< geom::cell >( "slab", 0,
    surfaces, std::vector<int>( { 1, -1 } ), std::make_shared<mc::material>( zr90 ) );

  const auto cells = std::vector< std::shared_ptr<geom::cell> >( { slab } );

  //create estimators
  auto leakage  = std::make_shared<mc::estimator>(
    "leakage",      0, mc::estimator_type_t::SURFACE_CURRENT,  std::nullopt, std::nullopt, std::nullopt, std::nullopt );
  std::vector< std::shared_ptr< mc::estimator > > estimators( { leakage } );

  //create source term
  std::function< mc::particle( util::RNG::ptr rng ) > src_dist = []( util::RNG::ptr rng ) {
    mc::particle p;
    p.pos = util::point( 1.0e-6, 0.0, 0.0 );
    p.dir = util::point( 1.0,    0.0, 0.0 );
    p.erg = 0.059;
    return p;
  };

  //create estmators for unresolved resonance sensitivity estimators
  auto rxn_ratios = mc::dosEstimator::rxn_ratio_vector(
    { std::make_pair( leakage, nullptr ) } );
  auto sens_params = std::vector< mc::dos_profile_t >(
     { mc::dos_profile_t( "zr90 total",   std::nullopt, 40090, xs::endf_mt_t::TOTAL,   std::nullopt ),
       mc::dos_profile_t( "zr90 capture", std::nullopt, 40090, xs::endf_mt_t::Z_GAMMA, std::nullopt ),
       mc::dos_profile_t( "zr90 elastic", std::nullopt, 40090, xs::endf_mt_t::ELASTIC, std::nullopt ) } );
  auto ur_sens = std::make_shared< mc::dosUREstimator >( xs, sens_params, rxn_ratios, mc::mc_calc_t::FIXED_SOURCE );
  const auto mc_calc = std::make_shared< mc::fixedSourceCalc >(
    num_histories, num_batches, src_dist, estimators, std::nullopt, ur_sens );

  //transport loop
  while ( ! mc_calc->is_done() ) {

    auto batch_size = mc_calc->batch_size();
    for ( auto history = 0 ; history < batch_size ; ++history ) {
       xs_buffer.flush();

       auto p = mc_calc->source_particle( history, rng );
       if ( ! p.cell ) {
         p.cell = geom::find_cell( p.pos, cells );
       }

       //variable to track number of collisions for analytical benchmark 3
       auto ncoll = 0;
       while ( p.alive() ) {
         //distance to collision
         const auto mat = p.cell->material();
         mat->fill_total_xs( p.erg, xs_buffer, rng );
         mat->fill_nufission_xs( p.erg, xs_buffer, rng );
         const auto dist_to_coll = xs_buffer.total_macro_xs > 0.0 ? -std::log( rng->sample() )/xs_buffer.total_macro_xs : std::numeric_limits<double>::max();

         //distance to boundary
         const auto next_surface = p.cell->next_surface( p.pos, p.dir );
         const auto dist_to_surf = next_surface.distance;

         const auto distance = std::fmin( dist_to_coll, dist_to_surf );
         p.move( distance );

         mc_calc->score( mc::estimator_type_t::TRACK_LENGTH_FLUX, p.wgt*distance, p.cell->id(), p, xs_buffer, rng );

         if ( distance == dist_to_coll ) {
           auto iso = mat->select_isotope( xs_buffer.total_macro_xs, xs_buffer.total_xs, rng );
           auto mt  = iso->select_rxn( p.erg, xs_buffer, rng );
           xs_buffer.collided_zaid = iso->zaid();
           xs_buffer.collided_mt   = mt;

           mc_calc->score( mc::estimator_type_t::COLLISION_FLUX, p.wgt/xs_buffer.total_macro_xs, p.cell->id(), p, xs_buffer, rng );

           if ( mt != xs::endf_mt_t::FISSION ) {
             auto yield = iso->yield( mt, p.erg );
             if ( yield > 0.0 ) {
               //auto [ Eout, dirout ] = iso->sample_secondary( mt, p.erg, p.dir, rng );
               //p.erg = Eout;
               //p.dir = dirout;
               //p.wgt *= yield;
               //special collision logic for analytical benchmark 3
               const auto r = rng->sample();
               p.erg = 0.054 + r * ( 0.059 - 0.054 );
               p.dir.x = r;
             }
             else {
               p.kill();
             }
           }
           else {
             //fission treated as capture
             p.kill();
           }
           ncoll++;
         }
         else {
           //hits surface
           if ( ncoll == 1 ) {
             mc_calc->score( mc::estimator_type_t::SURFACE_CURRENT, p.wgt, next_surface.surface->id(), p, xs_buffer, rng );
           }

           //check and process reflection if needed, then nudge by a small amount along surface normal
           //so the particle is in the next cell
           const auto& surf = next_surface.surface;
           if ( surf->bc() == geom::surface_bc_t::REFLECTING ) {
             p.dir = surf->reflect( p.pos, p.dir );
           }
           p.nudge( surf->normal( p.pos ) );
           p.cell = geom::find_cell( p.pos, cells );
           if ( ! p.cell ) p.kill();
         }
         //special logic for analytical benchmarks 3 to kill particle after two collisions
         if ( ncoll > 1 ) {
           p.kill();
         }
       }
       mc_calc->end_history();
    }
    mc_calc->end_batch();
  }

  //write analytical benchmark reference solutions
  xs[40090]->urr_data()->analytical_benchmark_3( 4.3675e-2, 10.0 );

  return 0;
}
