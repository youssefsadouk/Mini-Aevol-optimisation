//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <sstream>

Dna::Dna(int length, Threefry::Gen &&rng) {
    if (!shared_segment_) {
        // Create a new, independent sequence if no shared segment is being used
        seq_.resize(length);
        for (int32_t i = 0; i < length; i++) {
            seq_[i] = rng.random(NB_BASE) != 0; // Randomly generate DNA sequence
        }
    } else {
        // Use shared segment for DNA sequence
        assert(segment_end_ - segment_start_ == length); 
        // No need to populate seq_ as it will not be used
    }
}


int Dna::length() const {
    if (shared_segment_) {
        return segment_end_ - segment_start_;
    } else {
        return seq_.size();
    }
}

void Dna::save(gzFile backup_file) {
    std::stringstream bitset_stream;
    int dna_length;

    if (shared_segment_) {
        for (int i = segment_start_; i < segment_end_; ++i) {
            bitset_stream << shared_segment_->at(i);
        }
        dna_length = segment_end_ - segment_start_;
    } else {
        for (size_t i = 0; i < seq_.size(); ++i) {
            bitset_stream << seq_[i];
        }
        dna_length = seq_.size();
    }

    std::string bitset_string = bitset_stream.str();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, bitset_string.data(), dna_length * sizeof(char));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));
    std::string bitset_string(dna_length, '0');
    gzread(backup_file, &bitset_string[0], dna_length * sizeof(char));

    if (shared_segment_) {
        shared_segment_->lock();
        for (int i = 0; i < dna_length; ++i) {
            if (i + segment_start_ < shared_segment_->full_length()) {
                shared_segment_->set(i + segment_start_, bitset_string[i] == '1');
            } else {
                // Error handling: if the DNA data exceeds the bounds of the shared segment
                std::cerr << "Error: DNA data exceeds bounds of the shared segment." << std::endl;
                break;
            }
        }
        shared_segment_->unlock();
    } else {
        seq_.clear();
        for (char c : bitset_string) {
            seq_.push_back(c == '1');
        }
    }
}

bool Dna::is_part_of_shared_segment() const {
    return shared_segment_ != nullptr;
}

DnaSharedSegment* Dna::get_shared_segment() const {
    return shared_segment_;
}

void Dna::set_shared_segment(DnaSharedSegment *shared_segment, int start, int end) {
    shared_segment_ = shared_segment;
    segment_start_ = start;
    segment_end_ = end;
}

int Dna::get_segment_start() const {
    return segment_start_;
}

int Dna::get_segment_end() const {
    return segment_end_;
}

void Dna::set(int pos, char c) {
    if (shared_segment_) {
        if (pos >= segment_start_ && pos < segment_end_) {
            shared_segment_->set(pos - segment_start_, c == '1');
        }
    } else {
        assert(pos >= 0 && pos < seq_.size());
        seq_[pos] = (c == '1');
    }
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
    assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());
    if (shared_segment_) {
        if (pos_1 >= segment_start_ && pos_2 <= segment_end_) {
            // If the removal is entirely within the shared segment,
            // we need to update the shared segment and adjust the start and end of the current segment
            shared_segment_->lock();
            shared_segment_->remove(pos_1 - segment_start_, pos_2 - segment_start_);
            segment_end_ -= (pos_2 - pos_1);
            shared_segment_->unlock();
        }
    } else {
        int num_bits_to_remove = pos_2 - pos_1;
        for (int i = pos_2; i < seq_.size(); ++i) {
            seq_[i - num_bits_to_remove] = seq_[i];
        }
        seq_.resize(seq_.size() - num_bits_to_remove);
    }
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */

void Dna::insert(int pos, std::vector<char> seq) {
    assert(pos >= 0 && pos <= seq_.size());
    if (shared_segment_) {
        shared_segment_->lock();

        // Adjust position to relative position in shared segment
        int relative_pos = pos - segment_start_;
        assert(relative_pos >= 0 && relative_pos <= shared_segment_->full_length());

        // Perform insertion in the shared segment
        shared_segment_->insert(relative_pos, seq);

        shared_segment_->unlock();
    } else {
        seq_.resize(seq_.size() + seq.size());
        for (int i = seq_.size() - 1; i >= pos + seq.size(); --i) {
            seq_[i] = seq_[i - seq.size()];
        }
        for (size_t i = 0; i < seq.size(); ++i) {
            seq_[pos + i] = seq[i] == '1';
        }
    }
}


/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */

void Dna::insert(int pos, Dna* seq) {
    assert(pos >= 0 && pos <= seq_.size());

    if (shared_segment_) {
        shared_segment_->lock();

        // Adjust position to relative position in shared segment
        int relative_pos = pos - segment_start_;
        assert(relative_pos >= 0 && relative_pos <= shared_segment_->full_length());

        // Convert Dna segment to std::vector<char> for insertion
        std::vector<char> insert_seq(seq->seq_.size());
        for (size_t i = 0; i < seq->seq_.size(); ++i) {
            insert_seq[i] = seq->seq_[i] ? '1' : '0';
        }

        // Perform insertion in the shared segment
        shared_segment_->insert(relative_pos, insert_seq);

        shared_segment_->unlock();
    } else {
        // Original code for non-shared segment
        seq_.resize(seq_.size() + seq->seq_.size());
        for (int i = seq_.size() - 1; i >= pos + seq->seq_.size(); --i) {
            seq_[i] = seq_[i - seq->seq_.size()];
        }
        for (size_t i = 0; i < seq->seq_.size(); ++i) {
            seq_[pos + i] = seq->seq_[i];
        }
    }
}


void Dna::do_switch(int pos) {
    assert(pos >= 0 && pos < seq_.size());

    if (shared_segment_) {
        shared_segment_->lock();
        int relative_pos = pos - segment_start_;
        assert(relative_pos >= 0 && relative_pos < shared_segment_->full_length());

        // Perform the switch operation on the shared segment
        shared_segment_->flip(relative_pos);

        shared_segment_->unlock();
    } else {
        seq_.flip(pos);
    }
}


void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    assert(pos_1 >= 0 && pos_1 < seq_.size());
    assert(pos_2 >= 0 && pos_2 <= seq_.size());
    assert(pos_3 >= 0 && pos_3 <= seq_.size());

    if (shared_segment_) {
        shared_segment_->lock();

        // Adjust positions to relative positions in shared segment
        int relative_pos_1 = pos_1 - segment_start_;
        int relative_pos_2 = pos_2 - segment_start_;
        int relative_pos_3 = pos_3 - segment_start_;

        assert(relative_pos_1 >= 0 && relative_pos_1 < shared_segment_->full_length());
        assert(relative_pos_2 >= 0 && relative_pos_2 <= shared_segment_->full_length());
        assert(relative_pos_3 >= 0 && relative_pos_3 <= shared_segment_->full_length());

        // Perform the duplication operation on the shared segment
        shared_segment_->duplicate(relative_pos_1, relative_pos_2, relative_pos_3);

        shared_segment_->unlock();
    } else {
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
}


int Dna::promoter_at(int pos) {
    assert(pos >= 0 && (shared_segment_ ? pos < shared_segment_->full_length() : pos < seq_.size()));

    int dist_lead = 0;

    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        // Ensure we don't go past the end of the sequence
        int curr_pos = pos + motif_id;
        if (shared_segment_ ? curr_pos >= shared_segment_->full_length() : curr_pos >= seq_.size()) break;

        // Convert character in PROM_SEQ to a boolean
        bool expected_bit = PROM_SEQ[motif_id] == '1';

        // Access the bit either from shared_segment or from seq_
        bool actual_bit = shared_segment_ ? shared_segment_->at(curr_pos) == '1' : seq_[curr_pos];

        if (actual_bit != expected_bit) {
            dist_lead++;  
        }
    }

    return dist_lead;
}



// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
    assert(pos >= 0 && (shared_segment_ ? pos < shared_segment_->full_length() : pos < seq_.size()));

    int dist_term_lead = 0;

    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;
        int left = pos + (TERM_SIZE - 1) - motif_id;

        // If either 'right' or 'left' goes beyond the sequence length, stop checking
        if (shared_segment_ ? (right >= shared_segment_->full_length() || left >= shared_segment_->full_length()) : (right >= seq_.size() || left >= seq_.size())) break;

        // Check if the bits at 'right' and 'left' are different
        bool right_bit = shared_segment_ ? shared_segment_->at(right) == '1' : seq_[right];
        bool left_bit = shared_segment_ ? shared_segment_->at(left) == '1' : seq_[left];

        if (right_bit != left_bit) {
            dist_term_lead++; 
        }
    }

    return dist_term_lead;
}


bool Dna::shine_dal_start(int pos) {
    assert(pos >= 0 && (shared_segment_ ? pos < shared_segment_->full_length() : pos < seq_.size()));

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        int k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        int t_pos = pos + k_t;

        if (shared_segment_ ? t_pos >= shared_segment_->full_length() : t_pos >= seq_.size()) return false;
        
        char expected_char = SHINE_DAL_SEQ[k_t];

        if (expected_char != '*') {  
            bool expected_bit = expected_char == '1';

            bool actual_bit = shared_segment_ ? shared_segment_->at(t_pos) == '1' : seq_[t_pos];

            if (actual_bit != expected_bit) {
                return false;
            }
        }
    }

    return true;
}

bool Dna::protein_stop(int pos) {
    assert(pos >= 0 && (shared_segment_ ? pos < shared_segment_->full_length() : pos < seq_.size()));

    for (int k = 0; k < CODON_SIZE; k++) {
        int t_k = pos + k;

        if (shared_segment_ ? t_k >= shared_segment_->full_length() : t_k >= seq_.size()) {
            return false;
        }
        
        bool expected_bit = PROTEIN_END[k] == '1';

        bool actual_bit = shared_segment_ ? shared_segment_->at(t_k) == '1' : seq_[t_k];

        if (actual_bit != expected_bit) {
            return false;
        }
    }

    return true;
}

int Dna::codon_at(int pos) {
    assert(pos >= 0 && (shared_segment_ ? pos <= shared_segment_->full_length() - CODON_SIZE : pos <= seq_.size() - CODON_SIZE));

    int value = 0;

    for (int i = 0; i < CODON_SIZE; i++) {
        int t_pos = pos + i;

        if (shared_segment_ ? t_pos >= shared_segment_->full_length() : t_pos >= seq_.size()) {
            return -1;
        }

        bool bit = shared_segment_ ? shared_segment_->at(t_pos) == '1' : seq_[t_pos];

        if (bit) {
            value |= 1 << (CODON_SIZE - i - 1);
        }
    }

    return value;
}

