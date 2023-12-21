//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <sstream>

Dna::Dna(int length, Threefry::Gen &&rng) : seq_(length) {
    // Generate a random genome
    for (int32_t i = 0; i < length; i++) {
        seq_[i] = rng.random(NB_BASE) != 0; // Ensures that the result is a boolean
    }
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    // Create a stringstream for conversion
    std::stringstream bitset_stream;

    // Write the bitset to the stringstream
    bitset_stream << seq_;

    // Get the string from the stringstream
    std::string bitset_string = bitset_stream.str();

    // Calculate the length of the string
    int dna_length = bitset_string.length();

    // Write the length and the string to the file
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, bitset_string.data(), dna_length * sizeof(char));
}


void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    std::string bitset_string(dna_length, '0');
    gzread(backup_file, &bitset_string[0], dna_length * sizeof(char));

    seq_ = boost::dynamic_bitset<>(bitset_string);  // Construct bitset from string
}


void Dna::set(int pos, char c) {
    assert(pos >= 0 && pos < seq_.size());
    seq_[pos] = (c == '1');  // Sets the bit to true if c is '1', false otherwise
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
     assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());

     // Calculate the number of bits to be removed
     int num_bits_to_remove = pos_2 - pos_1;

     // Shift bits after pos_2 to the left
     for (int i = pos_2; i < seq_.size(); ++i) {
         seq_[i - num_bits_to_remove] = seq_[i];
     }

     // Resize the bitset to its new size
     seq_.resize(seq_.size() - num_bits_to_remove);
 }


/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */

void Dna::insert(int pos, std::vector<char> seq) {
    assert(pos >= 0 && pos <= seq_.size()); // Adjust to allow insertion at the end

    // Resize the bitset to accommodate the new bits
    seq_.resize(seq_.size() + seq.size());

    // Shift existing bits to the right to make room for the new sequence
    for (int i = seq_.size() - 1; i >= pos + seq.size(); --i) {
        seq_[i] = seq_[i - seq.size()];
    }

    // Insert the new bits
    for (size_t i = 0; i < seq.size(); ++i) {
        seq_[pos + i] = seq[i] == '1';
    }
}


/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */

void Dna::insert(int pos, Dna *seq) {
    assert(pos >= 0 && pos <= seq_.size()); // Adjust to allow insertion at the end

    // Resize the bitset to accommodate the new bits
    seq_.resize(seq_.size() + seq->seq_.size());

    // Shift existing bits to the right to make room for the new sequence
    for (int i = seq_.size() - 1; i >= pos + seq->seq_.size(); --i) {
        seq_[i] = seq_[i - seq->seq_.size()];
    }

    // Insert the new bits
    for (size_t i = 0; i < seq->seq_.size(); ++i) {
        seq_[pos + i] = seq->seq_[i];
    }
}

void Dna::do_switch(int pos) {
    assert(pos >= 0 && pos < seq_.size());
    seq_.flip(pos);  // Flip the bit at the specified position
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    assert(pos_1 >= 0 && pos_1 < seq_.size());
    assert(pos_2 >= 0 && pos_2 <= seq_.size());
    assert(pos_3 >= 0 && pos_3 <= seq_.size());

    boost::dynamic_bitset<> seq_dupl;
    size_t dupl_length = (pos_1 < pos_2) ? pos_2 - pos_1 : seq_.size() - pos_1 + pos_2;
    seq_dupl.resize(dupl_length);

    // Fill the duplicated segment bit by bit
    if (pos_1 < pos_2) {
        for (size_t i = 0; i < dupl_length; ++i) {
            seq_dupl[i] = seq_[pos_1 + i];
        }
    } else {
        size_t i = 0;
        // First part from pos_1 to end
        for (; pos_1 + i < seq_.size(); ++i) {
            seq_dupl[i] = seq_[pos_1 + i];
        }
        // Second part from start to pos_2
        for (size_t j = 0; j < pos_2; ++j, ++i) {
            seq_dupl[i] = seq_[j];
        }
    }

    // Resize the original bitset to accommodate the duplicated segment
    seq_.resize(seq_.size() + dupl_length);

    // Shift bits after pos_3 to make room for the duplicate segment
    for (int i = seq_.size() - dupl_length - 1; i >= pos_3; --i) {
        seq_[i + dupl_length] = seq_[i];
    }

    // Insert the duplicated segment
    for (size_t i = 0; i < dupl_length; ++i) {
        seq_[pos_3 + i] = seq_dupl[i];
    }
}


int Dna::promoter_at(int pos) {
    assert(pos >= 0 && pos < seq_.size());

    int dist_lead = 0;

    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        // Ensure we don't go past the end of the sequence
        if (pos + motif_id >= seq_.size()) break;

        // Convert character in PROM_SEQ to a boolean
        bool expected_bit = PROM_SEQ[motif_id] == '1';

        if (seq_[pos + motif_id] != expected_bit) {
            dist_lead++;  // Increment distance if there's a mismatch
        }
    }

    return dist_lead;
}


// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
    assert(pos >= 0 && pos < seq_.size());

    int dist_term_lead = 0;

    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;  
        int left = pos + (TERM_SIZE - 1) - motif_id;
         
		// If either 'right' or 'left' goes beyond the sequence length, stop checking
        if (right >= seq_.size() || left >= seq_.size()) break;

        // Check if the bits at 'right' and 'left' are different
        if (seq_[right] != seq_[left]) {
            dist_term_lead++;  // Increment if there's a mismatch
        }
    }

    return dist_term_lead;
}


bool Dna::shine_dal_start(int pos) {
    assert(pos >= 0 && pos < seq_.size());

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        int k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        int t_pos = pos + k_t;

		if (t_pos >= seq_.size()) return false;
		
        char expected_char = SHINE_DAL_SEQ[k_t];

        if (expected_char != '*') {  // Skip comparison for '*' characters
            bool expected_bit = expected_char == '1';

            if (seq_[t_pos] != expected_bit) {
                return false;
            }
        }
    }

    return true;
}


bool Dna::protein_stop(int pos) {
    assert(pos >= 0 && pos < seq_.size());

    for (int k = 0; k < CODON_SIZE; k++) {
        int t_k = pos + k;

		if (t_k >= seq_.size()) {
            return false;
        }
        
        // Convert character in PROTEIN_END to a boolean
        bool expected_bit = PROTEIN_END[k] == '1';

        if (seq_[t_k] != expected_bit) {
            return false;
        }
    }

    return true;
}


int Dna::codon_at(int pos) {
    assert(pos >= 0 && pos <= seq_.size() - CODON_SIZE);

    int value = 0;

    for (int i = 0; i < CODON_SIZE; i++) {
        int t_pos = pos + i;

        if (t_pos >= seq_.size()) {
            return -1;
        }

        if (seq_[t_pos]) {
            value |= 1 << (CODON_SIZE - i - 1);
        }
    }

    return value;
}
