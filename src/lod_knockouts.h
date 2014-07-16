
#ifndef _APOP_LOD_KNOCKOUTS_H_
#define _APOP_LOD_KNOCKOUTS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>

#include "gls.h"



namespace ealib {
    namespace analysis {
        template <typename EA>
        void setup_ko_ea(EA& parent_ea, EA& ea) {
            ea.rng().reset(get<RNG_SEED>(parent_ea));
            for(typename EA::population_type::iterator j=parent_ea.founder().begin(); j!=parent_ea.founder().end(); ++j) {
                
                typename EA::individual_ptr_type o1 = parent_ea.copy_individual(**j);
                o1->hw().initialize();
                ea.insert(ea.end(), o1);
            }
            
        }
        
        template <typename EA>
        void run_ea_res(EA& ea, int update_max) {
            int cur_update = 0;
            // and run till the group amasses the right amount of resources
            while ((get<GROUP_RESOURCE_UNITS>(ea,0) < get<GROUP_REP_THRESHOLD>(ea)) &&
                   (cur_update < update_max)){
                ea.update();
                ++cur_update;
            }
        }
        
        
        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         
         */
        LIBEA_ANALYSIS_TOOL(lod_knockouts_apops) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_knockouts.dat");
            df.add_field("lod_depth")
            .add_field("no_knockouts")
            .add_field("all_ifs_knockedout")
            .add_field("if_5_knockedout")
            .add_field("if_10_knockedout")
            .add_field("if_25_knockedout")
            .add_field("if_50_knockedout")
            .add_field("all_apops_knockedout")
            .add_field("apop_knockedout")
            .add_field("apop_5_knockedout")
            .add_field("apop_10_knockedout")
            .add_field("apop_25_knockedout")
            .add_field("apop_50_knockedout")
            .add_field("apop_x_knockedout");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control = ea.make_individual();
                setup_ko_ea(i->ea(), control->ea());
        
                typename EA::individual_ptr_type all_ifs_knockedout = ea.make_individual();
                // knockout all ifs
                knockout<if_workload_g5,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g10,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g25,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g50,instructions::nop_x>(all_ifs_knockedout->ea());
                setup_ko_ea(i->ea(), all_ifs_knockedout->ea());
                
                typename EA::individual_ptr_type if_5_knockedout = ea.make_individual();
                knockout<if_workload_g5,instructions::nop_x>(if_5_knockedout->ea());
                setup_ko_ea(i->ea(), if_5_knockedout->ea());

                typename EA::individual_ptr_type if_10_knockedout = ea.make_individual();
                knockout<if_workload_g10,instructions::nop_x>(if_10_knockedout->ea());
                setup_ko_ea(i->ea(), if_10_knockedout->ea());
                
                typename EA::individual_ptr_type if_25_knockedout = ea.make_individual();
                knockout<if_workload_g25,instructions::nop_x>(if_25_knockedout->ea());
                setup_ko_ea(i->ea(), if_25_knockedout->ea());
                
                typename EA::individual_ptr_type if_50_knockedout = ea.make_individual();
                knockout<if_workload_g50,instructions::nop_x>(if_50_knockedout->ea());
                setup_ko_ea(i->ea(), if_50_knockedout->ea());
                
                typename EA::individual_ptr_type all_apops_knockedout = ea.make_individual();
                knockout<instructions::apoptosis,instructions::nop_x>(all_apops_knockedout->ea());
                knockout<apop_g5,instructions::nop_x>(all_apops_knockedout->ea());
                knockout<apop_g10,instructions::nop_x>(all_apops_knockedout->ea());
                knockout<apop_g25,instructions::nop_x>(all_apops_knockedout->ea());
                knockout<apop_g50,instructions::nop_x>(all_apops_knockedout->ea());
                knockout<apop_gx,instructions::nop_x>(all_apops_knockedout->ea());
                setup_ko_ea(i->ea(), all_apops_knockedout->ea());
                
                typename EA::individual_ptr_type apop_knockedout = ea.make_individual();
                knockout<instructions::apoptosis,instructions::nop_x>(apop_knockedout->ea());
                setup_ko_ea(i->ea(), apop_knockedout->ea());
                
                typename EA::individual_ptr_type apop_5_knockedout = ea.make_individual();
                knockout<apop_g5,instructions::nop_x>(apop_5_knockedout->ea());
                setup_ko_ea(i->ea(), apop_5_knockedout->ea());
                
                typename EA::individual_ptr_type apop_10_knockedout = ea.make_individual();
                knockout<apop_g10,instructions::nop_x>(apop_10_knockedout->ea());
                setup_ko_ea(i->ea(), apop_10_knockedout->ea());
                
                typename EA::individual_ptr_type apop_25_knockedout = ea.make_individual();
                knockout<apop_g25,instructions::nop_x>(apop_25_knockedout->ea());
                setup_ko_ea(i->ea(), apop_25_knockedout->ea());
                
                typename EA::individual_ptr_type apop_50_knockedout = ea.make_individual();
                knockout<apop_g50,instructions::nop_x>(apop_50_knockedout->ea());
                setup_ko_ea(i->ea(), apop_50_knockedout->ea());

                typename EA::individual_ptr_type apop_x_knockedout = ea.make_individual();
                knockout<apop_gx,instructions::nop_x>(apop_x_knockedout->ea());
                setup_ko_ea(i->ea(), apop_x_knockedout->ea());
                
             
                run_ea_res(control->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(control->ea(), 0));
                
                run_ea_res(all_ifs_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(all_ifs_knockedout->ea(), 0));
                
                run_ea_res(if_5_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(if_5_knockedout->ea(), 0));
                
                run_ea_res(if_10_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(if_10_knockedout->ea(), 0));
                
                run_ea_res(if_25_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(if_25_knockedout->ea(), 0));
                
                run_ea_res(if_50_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(if_50_knockedout->ea(), 0));

                run_ea_res(all_apops_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(all_apops_knockedout->ea(), 0));
                
                run_ea_res(apop_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_knockedout->ea(), 0));
                
                run_ea_res(apop_5_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_5_knockedout->ea(), 0));

                run_ea_res(apop_10_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_10_knockedout->ea(), 0));
                
                run_ea_res(apop_25_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_25_knockedout->ea(), 0));
                
                run_ea_res(apop_50_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_50_knockedout->ea(), 0));
                
                run_ea_res(apop_x_knockedout->ea(), 1000);
                df.write(get<APOPTOSIS_COUNT>(apop_x_knockedout->ea(), 0));
                
                df.endl();
                
                
                ++lod_depth;
            }
            
            
        }

        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         
         */
        LIBEA_ANALYSIS_TOOL(lod_knockouts) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_knockouts.dat");
            df.add_field("lod_depth")
            .add_field("no_knockouts")
            .add_field("all_ifs_knockedout")
            .add_field("if_5_knockedout")
            .add_field("if_10_knockedout")
            .add_field("if_25_knockedout")
            .add_field("if_50_knockedout");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control = ea.make_individual();
                control->ea().rng().reset(get<RNG_SEED>(i->ea()));
                
                
                typename EA::individual_ptr_type all_ifs_knockedout = ea.make_individual();
                all_ifs_knockedout->ea().rng().reset(get<RNG_SEED>(i->ea()));
                // knockout all ifs
                knockout<if_workload_g5,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g10,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g25,instructions::nop_x>(all_ifs_knockedout->ea());
                knockout<if_workload_g50,instructions::nop_x>(all_ifs_knockedout->ea());
                
                
                typename EA::individual_ptr_type if_5_knockedout = ea.make_individual();
                if_5_knockedout->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<if_workload_g5,instructions::nop_x>(if_5_knockedout->ea());
                
                typename EA::individual_ptr_type if_10_knockedout = ea.make_individual();
                if_10_knockedout->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<if_workload_g10,instructions::nop_x>(if_10_knockedout->ea());
                
                typename EA::individual_ptr_type if_25_knockedout = ea.make_individual();
                if_25_knockedout->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<if_workload_g25,instructions::nop_x>(if_25_knockedout->ea());
                
                typename EA::individual_ptr_type if_50_knockedout = ea.make_individual();
                if_50_knockedout->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<if_workload_g50,instructions::nop_x>(if_50_knockedout->ea());
                

                
                // Setup founders!
                for(typename EA::individual_type::ea_type::population_type::iterator j=i->ea().founder().begin(); j!=i->ea().founder().end(); ++j) {
                
                    typename EA::individual_type::ea_type::individual_ptr_type o1 = i->ea().copy_individual(**j);
                    o1->hw().initialize();
                    control->ea().insert(control->ea().end(), o1);                    
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o2 = i->ea().copy_individual(**j);
                    o2->hw().initialize();
                    all_ifs_knockedout->ea().insert(all_ifs_knockedout->ea().end(), o2);
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o3 = i->ea().copy_individual(**j);
                    o3->hw().initialize();
                    if_5_knockedout->ea().insert(if_5_knockedout->ea().end(), o3);
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o4 = i->ea().copy_individual(**j);
                    o4->hw().initialize();
                    if_10_knockedout->ea().insert(if_10_knockedout->ea().end(), o4);
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o5 = i->ea().copy_individual(**j);
                    o5->hw().initialize();
                    if_25_knockedout->ea().insert(if_25_knockedout->ea().end(), o5);
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o6 = i->ea().copy_individual(**j);
                    o6->hw().initialize();
                    if_50_knockedout->ea().insert(if_50_knockedout->ea().end(), o6);
                }
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 1000;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control,0) < get<GROUP_REP_THRESHOLD>(*control)) &&
                       (cur_update < update_max)){
                    control->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(control->ea(), 0));
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*all_ifs_knockedout,0) < get<GROUP_REP_THRESHOLD>(*all_ifs_knockedout)) &&
                       (cur_update < update_max)){
                    all_ifs_knockedout->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(all_ifs_knockedout->ea(), 0));
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*if_5_knockedout,0) < get<GROUP_REP_THRESHOLD>(*if_5_knockedout)) &&
                       (cur_update < update_max)){
                    if_5_knockedout->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(if_5_knockedout->ea(), 0));
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*if_10_knockedout,0) < get<GROUP_REP_THRESHOLD>(*if_10_knockedout)) &&
                       (cur_update < update_max)){
                    if_10_knockedout->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(if_10_knockedout->ea(), 0));
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*if_25_knockedout,0) < get<GROUP_REP_THRESHOLD>(*if_25_knockedout)) &&
                       (cur_update < update_max)){
                    if_25_knockedout->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(if_25_knockedout->ea(), 0));
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*if_50_knockedout,0) < get<GROUP_REP_THRESHOLD>(*if_50_knockedout)) &&
                       (cur_update < update_max)){
                    if_50_knockedout->ea().update();
                    ++cur_update;
                }
                df.write(get<APOPTOSIS_COUNT>(if_50_knockedout->ea(), 0));
                
                
                df.endl();
                
                
                ++lod_depth;
            }
            
            
        }
        
        
    }
}
#endif
