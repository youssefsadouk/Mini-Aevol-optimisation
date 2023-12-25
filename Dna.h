//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <zlib.h>
#include "Threefry.h"
#include "aevol_constants.h"
#include <mutex>

class DnaSharedSegment;

class Dna {
public:
    Dna() = default;

    Dna(const Dna &clone) = default;

    Dna(int length, Threefry::Gen &&rng);

    ~Dna() = default;

    // Returns the length of the DNA segment
    int length() const;

    void save(gzFile backup_file);

    void load(gzFile backup_file);

    void set(int pos, char c);

    // Remove the DNA inbetween pos_1 and pos_2
    void remove(int pos_1, int pos_2);

    // Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, std::vector<char> seq);

    // Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, Dna *seq);

    void do_switch(int pos);

    void do_duplication(int pos_1, int pos_2, int pos_3);

    int promoter_at(int pos);

    int terminator_at(int pos);

    bool shine_dal_start(int pos);

    bool protein_stop(int pos);

    int codon_at(int pos);

    // Associate this DNA with a shared DNA segment
    void set_shared_segment(DnaSharedSegment *shared_segment, int start, int end);
    
	// Method to check if this DNA is part of a shared segment
    bool is_part_of_shared_segment() const;

    // Method to get the shared segment
    DnaSharedSegment* get_shared_segment() const;

    int get_segment_start() const;
    
    int get_segment_end() const;
    
    boost::dynamic_bitset<> seq_;

    // Pointer to shared DNA segment
    DnaSharedSegment *shared_segment_ = nullptr;
    int segment_start_ = 0;
    int segment_end_ = 0;
};

class DnaSharedSegment {
public:
    DnaSharedSegment(int total_length, Threefry::Gen &&rng) : full_sequence_(total_length) {
        // Initialize the full_sequence_ with random data
        for (int32_t i = 0; i < total_length; i++) {
            full_sequence_[i] = rng.random(NB_BASE) != 0;
        }
    }

    char at(int pos) const {
        return full_sequence_[pos] ? '1' : '0';
    }

    void set(int pos, bool value) {
        full_sequence_[pos] = value;
    }

    int full_length() const {
        return full_sequence_.size();
    }

    void lock() {
        // Lock the segment to prevent concurrent access
        mutex_.lock();
    }

    void unlock() {
        // Unlock the segment
        mutex_.unlock();
    }

    void remove(int pos_1, int pos_2) {
        assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= full_sequence_.size());
        int num_bits_to_remove = pos_2 - pos_1;
    
        // Create a temporary buffer to hold the updated sequence
        boost::dynamic_bitset<> temp_buffer(full_sequence_.size() - num_bits_to_remove);
    
        // Copy bits before the removal segment
        for (int i = 0; i < pos_1; ++i) {
            temp_buffer[i] = full_sequence_[i];
        }
    
        // Copy bits after the removal segment
        for (int i = pos_2; i < full_sequence_.size(); ++i) {
            temp_buffer[i - num_bits_to_remove] = full_sequence_[i];
        }
    
        // Replace the full sequence with the updated buffer
        full_sequence_ = std::move(temp_buffer);
    }

    void insert(int pos, const std::vector<char>& seq) {
        assert(pos >= 0 && pos <= full_sequence_.size());
        int seq_length = seq.size();
    
        // Create a temporary buffer to hold the updated sequence
        boost::dynamic_bitset<> temp_buffer(full_sequence_.size() + seq_length);
    
        // Copy bits before the insertion point
        for (int i = 0; i < pos; ++i) {
            temp_buffer[i] = full_sequence_[i];
        }
    
        // Insert the new sequence
        for (int i = 0; i < seq_length; ++i) {
            temp_buffer[pos + i] = seq[i] == '1';
        }
    
        // Copy the remaining bits after the insertion point
        for (size_t i = pos; i < full_sequence_.size(); ++i) {
            temp_buffer[i + seq_length] = full_sequence_[i];
        }
    
        // Replace the full sequence with the updated buffer
        full_sequence_ = std::move(temp_buffer);
    }

    void flip(int pos) {
        assert(pos >= 0 && pos < full_sequence_.size());
        full_sequence_.flip(pos);
    }

    void duplicate(int pos_1, int pos_2, int pos_3) {
        assert(pos_1 >= 0 && pos_1 < full_sequence_.size());
        assert(pos_2 >= 0 && pos_2 <= full_sequence_.size());
        assert(pos_3 >= 0 && pos_3 <= full_sequence_.size());
    
        size_t dupl_length = (pos_1 < pos_2) ? pos_2 - pos_1 : full_sequence_.size() - pos_1 + pos_2;
        
        boost::dynamic_bitset<> seq_dupl(dupl_length);
    
        // Fill the duplicated segment bit by bit
        if (pos_1 < pos_2) {
            for (size_t i = 0; i < dupl_length; ++i) {
                seq_dupl[i] = full_sequence_[pos_1 + i];
            }
        } else {
            size_t i = 0;
            for (; pos_1 + i < full_sequence_.size(); ++i) {
                seq_dupl[i] = full_sequence_[pos_1 + i];
            }
            for (size_t j = 0; j < pos_2; ++j, ++i) {
                seq_dupl[i] = full_sequence_[j];
            }
        }
    
        // Resize the original bitset to accommodate the duplicated segment
        full_sequence_.resize(full_sequence_.size() + dupl_length);
    
        // Shift bits after pos_3 to make room for the duplicate segment
        for (int i = full_sequence_.size() - dupl_length - 1; i >= pos_3; --i) {
            full_sequence_[i + dupl_length] = full_sequence_[i];
        }
    
        // Insert the duplicated segment
        for (size_t i = 0; i < dupl_length; ++i) {
            full_sequence_[pos_3 + i] = seq_dupl[i];
        }
    }
    
    
    boost::dynamic_bitset<> full_sequence_;
    mutable std::mutex mutex_;
};
