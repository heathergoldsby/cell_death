#ifndef _APOP_MOVIE_H_
#define _APOP_MOVIE_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>


namespace ealib {
    namespace analysis {
        
        
        /*! lod_movie
         */
        LIBEA_ANALYSIS_TOOL(movie) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::reverse_iterator i=lod.rbegin(); ++i;
            datafile df("movie.dat");
            typename EA::individual_ptr_type control_sp = ea.make_individual();
            control_sp->ea().rng().reset(get<RNG_SEED>((i->ea())));
            
            for(typename EA::individual_type::ea_type::population_type::iterator j=i->ea().founder().begin(); j!=i->ea().founder().end(); ++j) {
                typename EA::individual_type::ea_type::individual_ptr_type o = i->ea().copy_individual(**j);
                o->hw().initialize();
                control_sp->ea().insert(control_sp->ea().end(), o);
            }
            
            // hardcoded to run for 1000 updates.... can change
            int update_max = 1000;
            
            df.write(get<SPATIAL_X>(ea));
            df.write(get<SPATIAL_Y>(ea));
            df.endl();
            for (int i=0; i<=update_max; ++i) {
                control_sp->ea().update();
                df.write(i);
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                        typename EA::individual_type::ea_type::environment_type::location_ptr_type l = control_sp->ea().env().location(x,y);
                        
                        // Format: g(1)/s(0); fml; alive(1)/apop(0)
                        if (l->inhabitant() != 0) {
                            // l->occupied()

                            if (get<GERM_STATUS>(*(l->inhabitant()),true)) {
                                df.write("1");
                            } else {
                                df.write("0");
                            }
                            
                            df.write(get<WORKLOAD>(*(l->inhabitant()), 0));
                            
                            if (l->inhabitant()->alive()) {
                                df.write("1");
                            } else {
                                df.write("0");
                            }
                            
                            
                            
                        } else {
                            df.write("-1");
                            df.write("-1");
                            df.write("-1");
                        }
                        
                    }
                }
                df.endl();
                
            }
            
            df.endl();
            
        }
    }
}
#endif
