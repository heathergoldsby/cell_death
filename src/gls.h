//
//  gls.h
//  ealife
//
//  Created by Heather Goldsby on 8/23/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_GLS_H_
#define _EALIFE_GLS_H_


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>

#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>

#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/digital_evolution/extra_instruction_sets/cell_death.h>

#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ealib;
using namespace boost::accumulators;


/*
 Track results - tasks?
 Fixed inputs?
 
 */

LIBEA_MD_DECL(GERM_STATUS, "ea.gls.germ_status", int);
LIBEA_MD_DECL(STERILE, "ea.gls.sterile", int);
LIBEA_MD_DECL(REPLACEABLE, "ea.gls.replaceable", int);


LIBEA_MD_DECL(TASK_MUTATION_PER_SITE_P, "ea.gls.task_mutation_per_site_p", double);
//LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.gls.germ_mutation_per_site_p", double);
LIBEA_MD_DECL(WORKLOAD, "ea.gls.workload", double);
LIBEA_MD_DECL(TASK_MUTATION_MULT, "ea.gls.task_mutation_mult", double);
LIBEA_MD_DECL(NOT_MUTATION_MULT, "ea.gls.not_mutation_mult", double);
LIBEA_MD_DECL(NAND_MUTATION_MULT, "ea.gls.nand_mutation_mult", double);
LIBEA_MD_DECL(AND_MUTATION_MULT, "ea.gls.and_mutation_mult", double);
LIBEA_MD_DECL(ORNOT_MUTATION_MULT, "ea.gls.ornot_mutation_mult", double);
LIBEA_MD_DECL(OR_MUTATION_MULT, "ea.gls.or_mutation_mult", double);
LIBEA_MD_DECL(ANDNOT_MUTATION_MULT, "ea.gls.andnot_mutation_mult", double);
LIBEA_MD_DECL(NOR_MUTATION_MULT, "ea.gls.nor_mutation_mult", double);
LIBEA_MD_DECL(XOR_MUTATION_MULT, "ea.gls.xor_mutation_mult", double);
LIBEA_MD_DECL(EQUALS_MUTATION_MULT, "ea.gls.equals_mutation_mult", double);

LIBEA_MD_DECL(APOPTOSIS_COUNT, "ea.digevo.apoptosis_count", double);
LIBEA_MD_DECL(APOPTOSIS_WORKLOAD, "ea.digevo.apoptosis_workload", double);
LIBEA_MD_DECL(APOPTOSIS_SOMA_COUNT, "ea.digevo.apoptosis_soma_count", double);


// Germ instructions!

/*! Mark an organism as soma.
 */

DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
}

DIGEVO_INSTRUCTION_DECL(sterilize) {
    put<STERILE>(true,*p);
}

DIGEVO_INSTRUCTION_DECL(replaceable) {
    put<REPLACEABLE>(true,*p);
}

// checks to see if an organism is sterile before it divides
DIGEVO_INSTRUCTION_DECL(h_divide_respect_sterile) {
    if (!get<STERILE>(*p, false)) {
        if(hw.age() >= (0.8 * hw.original_size())) {
            typename Hardware::representation_type& r=hw.repr();
            
            // Check to see if the offspring would be a good length.
            int divide_pos = hw.getHeadLocation(Hardware::RH);
            int extra_lines = r.size() - hw.getHeadLocation(Hardware::WH);
            
            int child_size = r.size() - divide_pos - extra_lines;
            int parent_size = r.size() - child_size - extra_lines;
            double ratio = 2.0;
            
            if ((child_size < (hw.original_size()/ratio)) ||
                (child_size > (hw.original_size()*ratio)) ||
                (parent_size < (hw.original_size()/ratio)) ||
                (parent_size > (hw.original_size()*ratio))){
                // fail and die a miserable death!
                hw.replicated();
                return;
            }
            
            
            typename Hardware::representation_type::iterator f=r.begin(),l=r.begin();
            std::advance(f, hw.getHeadLocation(Hardware::RH));
            std::advance(l, hw.getHeadLocation(Hardware::WH));
            typename Hardware::representation_type offr(f, l);
            
            
            r.resize(parent_size);
            replicate(p, offr, ea);
            hw.replicated();
        }

    }
}

/*! Execute the next instruction if the organism is marked as germ.
 */

DIGEVO_INSTRUCTION_DECL(if_germ) {
    if(!get<GERM_STATUS>(*p,true)) {
        hw.advanceHead(Hardware::IP);
    }
}


/*! Execute the next instruction if the organism is marked as soma.
 */
DIGEVO_INSTRUCTION_DECL(if_soma){
    if(get<GERM_STATUS>(*p, true)) {
        hw.advanceHead(Hardware::IP);
    }
}

/*! Execute the next instruction if the wokload is greater than 5.
 */
DIGEVO_INSTRUCTION_DECL(if_workload_g5){
    if(get<WORKLOAD>(*p, 0) < 5) {
        hw.advanceHead(Hardware::IP);
    }
}

DIGEVO_INSTRUCTION_DECL(if_workload_g10){
    if(get<WORKLOAD>(*p, 0) < 10) {
        hw.advanceHead(Hardware::IP);
    }
}

DIGEVO_INSTRUCTION_DECL(if_workload_g25){
    if(get<WORKLOAD>(*p, 0) < 25) {
        hw.advanceHead(Hardware::IP);
    }
}

DIGEVO_INSTRUCTION_DECL(if_workload_g50){
    if(get<WORKLOAD>(*p, 0) < 50) {
        hw.advanceHead(Hardware::IP);
    }
}

/*! Execute the next instruction if the wokload is greater than 5.
 */
DIGEVO_INSTRUCTION_DECL(apop_g5){
    if(get<WORKLOAD>(*p, 0) > 5) {
        p->alive() = false;
        ea.events().death(*p,ea);
        put<APOPTOSIS_STATUS>(1, *p);
    }
}

DIGEVO_INSTRUCTION_DECL(apop_g10){
    if(get<WORKLOAD>(*p, 0) > 10) {
        p->alive() = false;
        ea.events().death(*p,ea);
        put<APOPTOSIS_STATUS>(1, *p);
    }
}

DIGEVO_INSTRUCTION_DECL(apop_g25){
    if(get<WORKLOAD>(*p, 0) > 25) {
        p->alive() = false;
        ea.events().death(*p,ea);
        put<APOPTOSIS_STATUS>(1, *p);
    }
}

DIGEVO_INSTRUCTION_DECL(apop_g50){
    if(get<WORKLOAD>(*p, 0) > 50) {
        p->alive() = false;
        ea.events().death(*p,ea);
        put<APOPTOSIS_STATUS>(1, *p);
    }
}

DIGEVO_INSTRUCTION_DECL(apop_gx){
    int rbx = hw.modifyRegister();

    if(get<WORKLOAD>(*p, 0) > hw.getRegValue(rbx)) {
        p->alive() = false;
        ea.events().death(*p,ea);
        put<APOPTOSIS_STATUS>(1, *p);
    }
}



/*! Selects the location of an empty neighbor location *or a neighbor marked as replacable!*
 to the parent as the location for an offspring. 
 (Note: here empty includes locations occupied by dead organisms.)
 
 If there is not an empty location, then the replacement does not proceed. This method
 does not makes sense with well-mixed, since the 'neighborhood' of an organism is
 8 random locations.
 */
struct empty_or_replaceable_neighbor {
    template <typename EA>
    std::pair<typename EA::environment_type::iterator, bool> operator()(typename EA::individual_ptr_type parent, EA& ea) {
        typedef typename EA::environment_type::iterator location_iterator;
        std::pair<location_iterator, location_iterator> i = ea.env().neighborhood(*parent);
        
        for( ; i.first != i.second; ++i.first) {
            if(!i.first->occupied() || get<REPLACEABLE>(*(i.first->inhabitant()),false)) {
                return std::make_pair(i.first, true);
            }
        }
        return std::make_pair(i.second, false);
    }
};


// Events!

/*! Did an organism die by suicide? If so, tally it!
 
 */

template <typename EA>
struct gs_apoptosis_event : death_event<EA> {
    
    //! Constructor.
    gs_apoptosis_event(EA& ea) : death_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~gs_apoptosis_event() {
    }
    
    /*! Called for every inheritance event. We are using the germ/soma status
     of the first parent
     */
    virtual void operator()(typename EA::individual_type& offspring,
                            EA& ea) {
        if (get<APOPTOSIS_STATUS>(offspring, 0) == 1) {
            get<APOPTOSIS_COUNT>(ea, 0) += 1;
            get<APOPTOSIS_WORKLOAD>(ea,0) += get<WORKLOAD>(offspring, 0);
            
            if (!get<GERM_STATUS>(offspring, true)) {
                get<APOPTOSIS_SOMA_COUNT>(ea,0) += 1;
            }
        }
        
    }
};

/*! An organism inherits its parent's germ/soma status. If it is undefined, 
 then it is set to germ.
 */
template <typename EA>
struct gs_inherit_event : inheritance_event<EA> {
    
    //! Constructor.
    gs_inherit_event(EA& ea) : inheritance_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~gs_inherit_event() {
    }
    
    /*! Called for every inheritance event. We are using the germ/soma status 
     of the first parent
     */
    virtual void operator()(typename EA::population_type& parents,
                            typename EA::individual_type& offspring,
                            EA& ea) {
        
        get<GERM_STATUS>(offspring, true) = get<GERM_STATUS>(**parents.begin(), true);
    }
};


/*! Triggers a task having a mutagenic effect on a organism.
 Configurable mutagenic rate for all tasks.
 */

template <typename EA>
struct task_mutagenesis : reaction_event<EA> {
    
    task_mutagenesis(EA& ea) : reaction_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::task_library_type::task_ptr_type task, // task pointer
                            double r, // amount of resource consumed
                            EA& ea) {

        double mult = get<TASK_MUTATION_MULT>(*task);
        double prob = get<TASK_MUTATION_PER_SITE_P>(ea) * mult;
        if (prob > 0) {
            configurable_per_site m(prob); 
            mutate(ind,m,ea);
            get<WORKLOAD>(ind,0.0) += mult;
        }
    }
};

/*! Triggers a task having a mutagenic effect on an organism in the colony with the same
 germ/soma status, but that has performed the least amount of work so far.
 Configurable mutagenic rate for all tasks.
 */

template <typename EA>
struct task_mutagenesis_control : reaction_event<EA> {
    
    task_mutagenesis_control(EA& ea) : reaction_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis_control() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::task_library_type::task_ptr_type task, // task pointer
                            double r, // amount of resource consumed
                            EA& ea) {
        
        
        double mult = get<TASK_MUTATION_MULT>(*task);
        double prob = get<TASK_MUTATION_PER_SITE_P>(ea) * mult;
        int gs_status = get<GERM_STATUS>(ind, true);
        typename EA::individual_type& sacrificial_org = ind;
        int smallest_workload = get<WORKLOAD>(ind, 0.0);
        if (prob > 0) {
            // search through the ea for an individual of the same g/s status, but
            // with the lowest workload
            for(typename EA::population_type::iterator j=ea.population().begin(); j!=ea.population().end(); ++j) {
                typename EA::individual_type& org=**j;
                if (get<GERM_STATUS>(org, true) == gs_status) {
                    if (get<WORKLOAD>(org,0.0) < smallest_workload) {
                        smallest_workload = get<WORKLOAD>(org,0.0);
                        sacrificial_org = org;
                    }
                }
            }
            
            configurable_per_site m(prob);
            mutate(sacrificial_org,m,ea);
            get<WORKLOAD>(sacrificial_org,0.0) += mult;
        }
    }
};

//! Performs group replication using germ lines.
template <typename EA>
struct gls_replication : end_of_update_event<EA> {
    //! Constructor.
    gls_replication(EA& ea) : end_of_update_event<EA>(ea), _df("gls.dat") {
        _df.add_field("update")
        .add_field("mean_germ_num")
        .add_field("mean_pop_num")
        .add_field("mean_germ_percent")
        .add_field("mean_germ_workload")
        .add_field("mean_germ_workload_var")
        .add_field("mean_soma_workload")
        .add_field("mean_soma_workload_var")
        .add_field("replication_count");
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~gls_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < 2) continue; 
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) && 
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                // grab a copy of the first individual: 
                typename EA::individual_type::ea_type::individual_type germ;
                int germ_present = false;
                
                // If so, setup a new replicate pop.
                // Find a germ...
                std::random_shuffle(i->population().begin(), i->population().end(), ea.rng());
                
                int germ_count = 0;
                int pop_count = 0;
                accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc; 
                accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc; 
                
                

                
                for(typename EA::individual_type::ea_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {

                    typename EA::individual_type::ea_type::individual_type& org=**j;
                    if (get<GERM_STATUS>(org, true)) {
                        ++germ_count;
                        germ_workload_acc(get<WORKLOAD>(org, 0.0));
                        if (!germ_present){
                            germ = org;
                            // Makes sure that we keep the size of the organism and discard its in-memory offspring
                            germ.repr().resize(org.hw().original_size());
                            germ.hw().initialize();
                            germ_present = true;

                        }
                    } else {
                        soma_workload_acc(get<WORKLOAD>(org, 0.0));
                    }
                    ++pop_count;
                }
                
                if (!germ_present) continue;
                
                pop_num.push_back(pop_count);
                germ_num.push_back(germ_count);
                germ_percent.push_back(germ_count/((double) i->population().size())*100.0); 
                germ_workload.push_back(mean(germ_workload_acc));
                germ_workload_var.push_back(variance(germ_workload_acc));
                
                if (germ_count != pop_count) {
                    soma_workload.push_back(mean(soma_workload_acc));
                    soma_workload_var.push_back(variance(soma_workload_acc));
                } else {
                    soma_workload.push_back(0);
                    soma_workload_var.push_back(0);
                }
                
                ++num_rep;
                
                
                if (germ_num.size() > 100) {
                    germ_num.pop_front();
                    germ_percent.pop_front();
                    pop_num.pop_front();
                    germ_workload.pop_front();
                    germ_workload_var.pop_front();
                    soma_workload.pop_front();
                    soma_workload_var.pop_front();
                }
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea)); 
                mutate(germ,m,p->ea());
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::ea_type::individual_ptr_type o=p->ea().copy_individual(germ.repr());
                inherits_from(germ, *o, p->ea());

                p->insert(p->end(), o);
                
                // add as founder
                //p->ea().founder().insert(p->ea().founder().end(), p->ea().copy_individual(*o));
                
                offspring.push_back(p);
                
                
                // reset resource units
                i->ea().env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*i.base());
                offspring_pop.push_back(p);
                inherits(parent_pop, offspring_pop, ea);
            }
        }
        
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<META_POPULATION_SIZE>(ea) - offspring.size(); 
            
            typename EA::population_type survivors;
            select_n<selection::random< > >(ea.population(), survivors, n, ea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(ea.population(), survivors);
        }
        
        
        if ((ea.current_update() % 100) == 0) {
            if (germ_num.size() > 0) {
                _df.write(ea.current_update())
                .write(std::accumulate(germ_num.begin(), germ_num.end(), 0.0)/germ_num.size())
                .write(std::accumulate(pop_num.begin(), pop_num.end(), 0.0)/pop_num.size())
                .write(std::accumulate(germ_percent.begin(), germ_percent.end(), 0.0)/germ_percent.size())
                .write(std::accumulate(germ_workload.begin(), germ_workload.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(germ_workload_var.begin(), germ_workload_var.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(soma_workload.begin(), soma_workload.end(), 0.0)/soma_workload.size())
                .write(std::accumulate(soma_workload_var.begin(), soma_workload_var.end(), 0.0)/soma_workload.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
            } else {
                _df.write(ea.current_update())
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .endl();     
            }
        }
    }
    
    datafile _df;    
    std::deque<double> germ_num; 
    std::deque<double> germ_percent;
    std::deque<double> pop_num;
    std::deque<double> germ_workload;
    std::deque<double> germ_workload_var;
    std::deque<double> soma_workload;
    std::deque<double> soma_workload_var;
    int num_rep;
    
    
};

/*! Prints information about apoptotic cells
 */


template <typename EA>
struct apoptosis_tracking : end_of_update_event<EA> {
    apoptosis_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("apop.dat") {
        _df.add_field("update")
        .add_field("mean_apop")
        .add_field("max_apop")
        .add_field("mean_apop_workload")
        .add_field("max_apop_workload")
        .add_field("mean_apop_soma")
        .add_field("max_apop_soma");
        
    }
    
    //! Destructor.
    virtual ~apoptosis_tracking() {
    }
    
    //! Track how many task-switches are being performed!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            accumulator_set<double, stats<tag::mean, tag::max> > apop;
            accumulator_set<double, stats<tag::mean, tag::max> > apop_wl;
            accumulator_set<double, stats<tag::mean, tag::max> > apop_s;

            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                apop(get<APOPTOSIS_COUNT>(i->ea(), 0));
                apop_wl(get<APOPTOSIS_WORKLOAD>(i->ea(), 0));
                apop_s(get<APOPTOSIS_SOMA_COUNT>(i->ea(),0));
            }

            _df.write(ea.current_update())
            .write(mean(apop))
            .write(max(apop))
            .write(mean(apop_wl))
            .write(max(apop_wl))
            .write(mean(apop_s))
            .write(max(apop_s))
            .endl();
    
        }
        
    }
    datafile _df;
    
};



#endif
