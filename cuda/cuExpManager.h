//
// Created by elturpin on 03/12/2020.
//

#pragma once

#include "Abstract_ExpManager.h"
#include "ExpManager.h"
#include "RandService.h"

class cuIndividual;

class cuExpManager : public Abstract_ExpManager {
public:
    explicit cuExpManager(const ExpManager* cpu_exp);

    ~cuExpManager() override;

    void run_evolution(int nb_gen) override;

    void save(int t) const final;

    void load(int t) final;

private:
    void run_a_step();
    void evaluate_population();

    // Interface Host - Device
    void transfer_to_device();
    void transfer_to_host() const;
    void device_data_destructor();

    // Host Data
    int nb_indivs_;

    int genome_length_;
    char** host_individuals_;  // Consider replacing with a single large DNA structure

    key_value_type seed_;
    size_t nb_counter_;
    ctr_value_type* counters_;

    double* target_;

    int grid_height_;
    int grid_width_;

    double mutation_rate_;

    int backup_step_;

    // Device data
    cuIndividual* device_individuals_;
    char* all_child_genome_;  // Consider using segments of a larger DNA structure
    char* all_parent_genome_;  // Consider using segments of a larger DNA structure

    int* reproducers_;

    double* device_target_;

    RandService* rand_service_;

    // DNA Segment Information
    uint dna_segment_start_;  // Start position of the segment in the large DNA structure
    uint dna_segment_length_;  // Length of the segment in the large DNA structure
};


