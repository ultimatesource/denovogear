/*
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#ifndef CXX_HTS_BAM_H
#define CXX_HTS_BAM_H

#include "hts.h"

#include <htslib/sam.h>

namespace hts {
namespace bam {

typedef bam1_t BareAlignment;

class File;

typedef std::pair<const uint32_t *, const uint32_t *> cigar_t;
typedef std::pair<const uint8_t *, const uint8_t *> data_t;

class Alignment : protected BareAlignment {
public:
    Alignment() {
        data = nullptr;
        m_data = 0;
        l_data = 0;
    }

    Alignment(const Alignment &other) {
        data = nullptr;
        m_data = 0;
        bam_copy1(base(), other.base());
    }

    Alignment(Alignment &&other) : BareAlignment(*other.base()) {
        other.data = nullptr;
        other.l_data = other.m_data = 0;
    }

    ~Alignment() {
        free(data);
    }

    Alignment &operator=(const Alignment &other) {
        if(this != &other) {
            bam_copy1(base(), other.base());
        }
        return *this;
    }

    Alignment &operator=(Alignment &&other) {
        if(this == &other) {
            return *this;
        }
        free(data);
        *base() = *other.base();
        other.data = nullptr;
        other.l_data = other.m_data = 0;
        return *this;
    }

    inline uint32_t target_id() const { return core.tid; }
    inline uint32_t position() const { return core.pos; }
    inline uint32_t map_qual() const { return core.qual; }
    inline uint32_t mate_target_id() const { return core.mtid; }
    inline uint32_t mate_position() const { return core.mpos; }
    inline uint16_t flags() const { return core.flag; }

    inline bool is_any(uint16_t f) const { return ((flags() & f) != 0); }
    inline bool are_only(uint16_t f) const { return (flags() == f); }
    inline bool are_all(uint16_t f) const { return ((flags() & f) == f); }

    inline bool is_reversed() const { return is_any(BAM_FREVERSE); }
    inline bool mate_is_reversed() const { return is_any(BAM_FMREVERSE); }

    inline const char *qname() const {
        return bam_get_qname(this);
    }
    inline cigar_t cigar() const {
        const uint32_t *b = bam_get_cigar(this);
        return std::make_pair(b, b + core.n_cigar);
    }
    //TODO: Make seq-specific iterator
    inline data_t seq() const {
        const uint8_t *b = bam_get_seq(this);
        return std::make_pair(b, b + (core.l_qseq + 1) / 2);
    }
    inline data_t seq_qual() const {
        const uint8_t *b = bam_get_qual(this);
        return std::make_pair(b, b + core.l_qseq);
    }
    inline data_t aux() const {
        const uint8_t *b = bam_get_aux(this);
        return std::make_pair(b, b + bam_get_l_aux(this));
    }

    inline uint8_t seq_at(int32_t x) const {
        assert(0 <= x && x < core.l_qseq);
        return bam_seqi(bam_get_seq(this), x);
    }

    inline uint8_t *aux_get(const char tag[3]) {
        return bam_aux_get(base(), tag);
    }

protected:
    BareAlignment *base() {return static_cast<BareAlignment *>(this);}
    const BareAlignment *base() const {return static_cast<const BareAlignment *>(this);}

    friend class File;
};

class File : public hts::File {
public:
    File(hts::File &&other, const char *region = nullptr,
         const char *fasta = nullptr,
         int min_mapQ = 0) : hts::File(std::move(other)),
        hdr_{nullptr, bam_hdr_destroy}, iter_{nullptr, hts_itr_destroy},
        min_mapQ_(min_mapQ) {
        // TODO: handle different modes here???
        // TODO: remove throws or move some of them to base class???

        if(!is_open()) { // Return early if the opening failed.  User must also check is_open to detect failure.
            return;
        }
        if(format().category != sequence_data)
            throw std::runtime_error("file '" + std::string(name())
                                     + "' does not contain sequence data (BAM/SAM/CRAM).");

        SetFaiFileName(fasta);
        hdr_.reset(sam_hdr_read(handle()));
        if(!hdr_) {
            throw std::runtime_error("unable to read header in file '" + std::string(
                                         name()) + "'.");
        }
        if(region != nullptr && region[0] != '\0') {
            std::unique_ptr<hts_idx_t, void(*)(hts_idx_t *)> idx{
                sam_index_load(handle(), name()), hts_idx_destroy};
            if(!idx) {
                throw std::runtime_error("unable to load index for '" + std::string(
                                             name()) + "'.");
            }
            iter_.reset(sam_itr_querys(idx.get(), hdr_.get(), region));

            if(!iter_) {
                throw std::runtime_error("unable to parse region '" + std::string(region) +
                                         "' in file '" + std::string(name()) + "'.");
            }
        }
    }

    File(const char *file, const char *mode, const char *region = nullptr,
         const char *fasta = nullptr,
         int min_mapQ = 0) : File(hts::File(file, mode), region, fasta, min_mapQ) {
    }

    int Read(Alignment *p) {
        int ret;
        for(;;) {
            ret = (iter_) ? sam_itr_next(handle(), iter_.get(), p->base())
                  : sam_read1(handle(), hdr_.get(), p->base());
            if(ret < 0) {
                break;
            }
            if(p->is_any(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
                continue;
            }
            if(p->map_qual() < min_mapQ_) {
                continue;
            }
            break;
        }
        return ret;
    }
    int operator()(Alignment *p) { return Read(p); }

    int Write(const bam1_t &b) {
        return sam_write1(handle(), hdr_.get(), &b);
    }

    const bam_hdr_t *header() const { return hdr_.get(); }
    const hts_itr_t *iter() const { return iter_.get(); }

protected:
    std::unique_ptr<bam_hdr_t, void(*)(bam_hdr_t *)> hdr_; // the file header
    std::unique_ptr<hts_itr_t, void(*)(hts_itr_t *)>
    iter_; // NULL if a region not specified

    int min_mapQ_; // mapQ filter
};

} // namespace bam

} // namespace hts

#endif