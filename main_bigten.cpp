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
  //eigenvalue simulation parameters
  const auto batch_size_initial = 50000;
  const auto num_inactive =  20;
  const auto num_cycles   =  270;

  //cross sections and materials
  std::vector< std::shared_ptr<xs::aceData> > ace;


  // big-ten data
  ace.push_back( std::make_shared< xs::aceData >( "data/92234.710nc" ) );
  ace.push_back( std::make_shared< xs::aceData >( "data/92235.710nc" ) );
  ace.push_back( std::make_shared< xs::aceData >( "data/92236.710nc" ) );
  ace.push_back( std::make_shared< xs::aceData >( "data/92238.710nc" ) );

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

  //big-ten, keff = xxxxxxx
  auto ieu10 = std::make_shared<mc::material>();
  ieu10->add( 4.8416e-05, xs[92234] );
  ieu10->add( 4.8151e-03, xs[92235] );
  ieu10->add( 1.7407e-05, xs[92236] );
  ieu10->add( 4.3181e-02, xs[92238] );

  auto depleted_u = std::make_shared<mc::material>();
  depleted_u->add( 2.8672e-07, xs[92234] );
  depleted_u->add( 1.0058e-04, xs[92235] );
  depleted_u->add( 1.1468e-06, xs[92236] );
  depleted_u->add( 4.7677e-02, xs[92238] );

  const auto inner_cyl = std::make_shared< geom::cylinder_z >( "inner cyl", 0, 0.0, 0.0, 26.67 );
  const auto outer_cyl = std::make_shared< geom::cylinder_z >( "outer cyl", 1, 0.0, 0.0, 41.91 );
  const auto zplane1   = std::make_shared< geom::plane >( "zplane1", 2, 0.0, 0.0, 1.0, -48.2600 );
  const auto zplane2   = std::make_shared< geom::plane >( "zplane2", 3, 0.0, 0.0, 1.0, -28.8169 );
  const auto zplane3   = std::make_shared< geom::plane >( "zplane3", 4, 0.0, 0.0, 1.0, +28.8169 );
  const auto zplane4   = std::make_shared< geom::plane >( "zplane4", 5, 0.0, 0.0, 1.0, +48.2600 );

  const auto core_surfaces = std::vector< std::shared_ptr<geom::surface> >(
   { inner_cyl, zplane2, zplane3 } );
  const auto reflector_bottom_surfaces = std::vector< std::shared_ptr<geom::surface> >(
   { outer_cyl, zplane1, zplane2 } );
  const auto reflector_middle_surfaces = std::vector< std::shared_ptr<geom::surface> >(
   { inner_cyl, outer_cyl, zplane2, zplane3 } );
  const auto reflector_top_surfaces = std::vector< std::shared_ptr<geom::surface> >(
   { outer_cyl, zplane3, zplane4 } );

  const auto core             = std::make_shared< geom::cell >(
    "core",             0, core_surfaces,             std::vector<int>( {-1, +1, -1}     ), ieu10 );
  const auto reflector_bottom = std::make_shared< geom::cell >(
    "bottom reflector", 1, reflector_bottom_surfaces, std::vector<int>( {-1, +1, -1}     ), depleted_u );
  const auto reflector_middle = std::make_shared< geom::cell >(
    "middle reflector", 2, reflector_middle_surfaces, std::vector<int>( {+1, -1, +1, -1} ), depleted_u );
  const auto reflector_top    = std::make_shared< geom::cell >(
    "top reflector",    3, reflector_top_surfaces,    std::vector<int>( {-1, +1, -1}     ), depleted_u );

  const auto cells = std::vector< std::shared_ptr<geom::cell> >( { core, reflector_bottom, reflector_middle, reflector_top } );

  auto leakage  = std::make_shared<mc::estimator>(
    "leakage", 0, mc::estimator_type_t::SURFACE_CURRENT, std::set<uint32_t>({ 1, 2, 5}), std::nullopt, std::nullopt, std::nullopt );
  auto keff     = std::make_shared<mc::estimator>(
    "keff",    1, mc::estimator_type_t::TRACK_LENGTH_FLUX, std::nullopt, std::nullopt, xs::endf_mt_t::NUFISSION, std::nullopt );
  auto core_u238_capture = std::make_shared<mc::estimator>(
    "core u238 capture", 2, mc::estimator_type_t::TRACK_LENGTH_FLUX, std::set<uint32_t>( {0} ), 92238, xs::endf_mt_t::Z_GAMMA, std::nullopt );
  auto core_u235_fission = std::make_shared<mc::estimator>(
    "core u235 fission", 3, mc::estimator_type_t::TRACK_LENGTH_FLUX, std::set<uint32_t>( {0} ), 92235, xs::endf_mt_t::FISSION, std::nullopt );
  std::vector< std::shared_ptr< mc::estimator > > estimators( { leakage, keff, core_u238_capture, core_u235_fission } );

  const auto ebinning = std::vector<double>( { 0.0, 2.0000e-02, 1.4900e-01, 30.0 } );
  auto rxn_ratios = mc::dosEstimator::rxn_ratio_vector(
    { std::make_pair( leakage, keff ), std::make_pair( core_u238_capture, core_u235_fission ) } );
  auto sens_params = std::vector< mc::dos_profile_t >(
     { mc::dos_profile_t( "u235 fission", std::nullopt, 92235, xs::endf_mt_t::FISSION, std::nullopt ),
       mc::dos_profile_t( "u238 capture", std::nullopt, 92238, xs::endf_mt_t::Z_GAMMA, std::nullopt ) } );
  auto ur_sens = std::make_shared< mc::dosUREstimator >( xs, sens_params, rxn_ratios, mc::mc_calc_t::K_EIGENVALUE );

  const auto mc_calc = std::make_shared< mc::kEigenvalueCalc >(
    batch_size_initial, num_inactive, num_cycles, util::point( 0.0, 0.0, 0.0 ), cells, estimators, std::nullopt, ur_sens );

  while ( ! mc_calc->is_done() ) {

    auto batch_size = mc_calc->batch_size();
    for ( auto history = 0 ; history < batch_size ; ++history ) {
       xs_buffer.flush();

       auto p = mc_calc->source_particle( history, rng );
       if ( ! p.cell ) {
         p.cell = geom::find_cell( p.pos, cells );
       }
       //keff_sens.start_history( source_bank[history].dos_id );

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
               auto [ Eout, dirout ] = iso->sample_secondary( mt, p.erg, p.dir, rng );
               p.erg = Eout;
               p.dir = dirout;
               p.wgt *= yield;
             }
             else {
               p.kill();
             }
           }
           else {
             //fission treated as capture
             p.kill();
           }
         }
         else {
           //hits surface
           mc_calc->score( mc::estimator_type_t::SURFACE_CURRENT, p.wgt, next_surface.surface->id(), p, xs_buffer, rng );

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
       }
       mc_calc->end_history();
    }
    mc_calc->end_batch();
  }
  return 0;
}
