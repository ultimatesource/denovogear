/*
 * Copyright (c) 2015 Reed A. Cartwright <reed@cartwrig.ht>
 * Copyright (c) 2015 Kael Dai <kdai1@asu.edu>
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
#ifndef CXX_HTS_BCF_H
#define CXX_HTS_BCF_H

#include "hts.h"

#include <vector>
#include <string>
#include <set>
#include <htslib/vcf.h>

extern "C" {
    // The htslib header does not match the library
    void bcf_empty1(bcf1_t *v);
}

namespace hts {
namespace bcf {

const float float_missing = []() -> float {
    // We use a lambda function so the result can be constant.
    float a;
    bcf_float_set_missing(a);
    return a;
}();
const int32_t int32_missing = bcf_int32_missing;
const int16_t int16_missing = bcf_int16_missing;
const int8_t int8_missing = bcf_int8_missing;

typedef bcf1_t BareVariant;

class File;

class Variant : protected BareVariant {
public:
    Variant() = default;

    Variant(const File &file);

    ~Variant() {
        bcf_empty1(base());
    }

    void Clear() {
        bcf_clear(base());
    }

    // Getters
    int32_t target_id() const { return rid; }
    int32_t position() const { return pos; }
    int32_t ref_length() const { return rlen; }
    float quality() const { return qual; }

    // Setters
    void quality(float q) { qual = q; }
    void target(const char *chrom) {
        rid = (chrom != nullptr) ? bcf_hdr_name2id(hdr(), chrom) : -1;
    }
    void position(int32_t p) { pos = p; }
    void id(const char *str) {
        bcf_update_id(hdr(), base(), str);
    }
    int filter(const char *str) {
        if(str == nullptr) {
            return 0;
        }
        int32_t fid = bcf_hdr_id2int(hdr(), BCF_DT_ID, str);
        return bcf_add_filter(hdr(), base(), fid);
    }
    int filter(const std::string &str) {
        int32_t fid = bcf_hdr_id2int(hdr(), BCF_DT_ID, str.c_str());
        return bcf_add_filter(hdr(), base(), fid);
    }

    /**
     * alleles() - Set the REF and ALT fields in the current record
     * @str: A comma separated list of all the alleles that show up in the sample.
     *       the first element is the REF value.
     */
    int alleles(const std::string &str) {
        return bcf_update_alleles_str(hdr(), base(), str.c_str());
    }
    int alleles(const char *str) {
        if(str == nullptr) {
            return 0;
        }
        return bcf_update_alleles_str(hdr(), base(), str);
    }

    /**
     * info() - Add another key=Value pair to the INFO field
     * @key: must be defined in the Header or won't be added.
     * @value: if set to NULL then previous set key-value pairs will be removed.
     */
    int info(const char *key, const float *value, std::size_t count) {
        assert(key != nullptr);
        return bcf_update_info_float(hdr(), base(), key, value, count);
    }
    int info(const char *key, const int32_t *value, std::size_t count) {
        assert(key != nullptr);
        return bcf_update_info_int32(hdr(), base(), key, value, count);
    }
    int info(const char *key, const std::string &value) {
        assert(key != nullptr);
        return bcf_update_info_string(hdr(), base(), key, value.c_str());
    }
    int info(const char *key, const char *value) {
        assert(key != nullptr);
        return bcf_update_info_string(hdr(), base(), key, value);
    }

    template<typename T>
    int info(const char *key, T value) {
        return info(key, &value, 1);
    }
    template<typename T, typename A>
    int info(const char *key, const std::vector<T, A> &value) {
        return info(key, &value[0], value.size());
    }
    template<typename T, std::size_t N>
    int info(const char *key, const std::array<T, N> &value) {
        return info(key, &value[0], value.size());
    }
    template<typename T, std::size_t N>
    int info(const char *key, const T(&value)[N]) {
        return info(key, &value[0], N);
    }

    /**
     * samples() - Add a key to the FORMAT field and update each sample column
     *                  with the corresponding key values.
     * @name: the tag name that appears in FORMAT field. Must be defined in the header
     * @data: A list of values that will populate the sample/genotype fields. The vector
     *        should be a multiple of the number sample columns.
     */
    int samples(const char *name, const std::vector<float> &data) {
        assert(name != nullptr);
        return bcf_update_format_float(hdr(), base(), name, &data[0], data.size());
    }
    int samples(const char *name, const std::vector<int32_t> &data) {
        assert(name != nullptr);
        return bcf_update_format_int32(hdr(), base(), name, &data[0], data.size());
    }
    int samples(const char *name, const std::vector<const char *> &data) {
        assert(name != nullptr);
        return bcf_update_format_string(hdr(), base(), name,
                                        const_cast<const char **>(&data[0]), data.size());
    }
    int samples(const char *name, const std::vector<std::string> &data) {
        std::vector<const char *> v(data.size());
        for(decltype(data.size()) u = 0; u < data.size(); ++u) {
            v[u] = data[u].c_str();
        }
        return samples(name, v);
    }
    int sample_genotypes(const std::vector<int32_t> &data) {
        return bcf_update_genotypes(hdr(), base(), &data[0], data.size());
    }

protected:
    BareVariant *base() {return static_cast<BareVariant *>(this);}
    const BareVariant *base() const {return static_cast<const BareVariant *>(this);}

    const bcf_hdr_t *hdr() {
        // check if the header pointer has not been initialized
        assert(hdr_);
        return hdr_.get();
    }

private:
    std::shared_ptr<const bcf_hdr_t> hdr_;

    friend class File;
};

/**
 * File - Class for writing VCF/BCF files to the output, either stdout or to a file
 * To write to stdout set file="-" in constructor. To write as a BCF file set mode="wb" in constructor.
 * To use:
 *  1. call AddHeaderMetadata() and AddSample() to add information to the VCF header
 *  2. call WriteHdr() to finish the header
 *  3. populate the first few columns of a VCF body line using SetID(), SetFilter(), SetAlleles(), UpdateInfo() and UpdateSamples()
 *  4. call WriteRecord() to save/print the current VCF line and move on to the next
 */
class File : public hts::File {
public:

    /**
     * Initialize File/output stream for writing VCF/BCF
     * @file - file name to write to, use "-" for stdout
     * @mode: "w" = write VCF, "wb" = write BCF
     * @source: name of application writing the VCF file (optional)
     */
    File(hts::File &&other) : hts::File(std::move(other)) {
        if(!is_open()) { // nothing to do
            return;
        }

        if(is_write()) {
            hdr_ = std::shared_ptr<bcf_hdr_t>(bcf_hdr_init("w"), bcf_hdr_destroy);
        } else {
            if(format().category != variant_data)
                throw std::runtime_error("file '" + std::string(name())
                                         + "' does not contain variant data (BCF/VCF).");
            hdr_ = std::shared_ptr<bcf_hdr_t>(bcf_hdr_read(handle()), bcf_hdr_destroy);
        }
        if(!hdr_) {
            throw std::runtime_error("unable to access header in file '" + std::string(
                                         name()) + "'.");
        }
    }

    File(const char *file, const char *mode) : File(hts::File(file, mode)) {
    }

    Variant InitVariant() const {
        return Variant{*this};
    }

    /** Use this version to add a new INFO/FILTER/FORMAT string to the header */
    int AddHeaderMetadata(const char *line) {
        if(line == nullptr) {
            return -1;    //bcf_hdr_append returns if line cannot be processed
        }
        return bcf_hdr_append(hdr(), line);
    }
    int AddHeaderMetadata(const std::string &line) {
        return bcf_hdr_append(hdr(), line.c_str());
    }

    /** AddHeaderMetadata() - Adds a "##key=value" line to the VCF header */
    int AddHeaderMetadata(const char *key, const char *value) {
        if(key == nullptr || value == nullptr) {
            return -1;
        }
        std::string line = std::string("##") + key + "=" + value;
        return bcf_hdr_append(hdr(), line.c_str());
    }

    int AddHeaderMetadata(const char *key, const std::string &value) {
        return AddHeaderMetadata(key, value.c_str());
    }

    template<typename T>
    int AddHeaderMetadata(const char *key, T value) {
        return AddHeaderMetadata(key, std::to_string(value));
    }

    /** Add another sample/genotype field to the VCF file. */
    int AddSample(const char *sample) {
        return bcf_hdr_add_sample(hdr(), sample);
    }

    /** Add a "#contig=" metadata entry to the VCF header. */
    int AddContig(const char *contigid, uint32_t length) {
        if(contigid == nullptr) {
            return -1;
        }
        std::string conv = std::string("##contig=<ID=") + contigid
                           + ",length=" + std::to_string(length) + ">";
        return bcf_hdr_append(hdr(), conv.c_str());
    }

    /** Creates the header field. Call only after adding all the sample fields */
    int WriteHeader() {
        // htslib requires NULL sample before it will write out all the other samples
        bcf_hdr_add_sample(hdr(), nullptr);
        return bcf_hdr_write(handle(), hdr());
    }

    std::pair<char **, int> samples() const {
        return {hdr_->samples, bcf_hdr_nsamples(header())};
    }

    const bcf_hdr_t *header() const { return hdr_.get(); }

    /** Writes out the up-to-date info in the record and prepares for the next line */
    void WriteRecord(Variant &rec) {
        // Add line to the body of the VCF
        assert(rec.hdr() == hdr());
        bcf_write(handle(), hdr(), &rec);
    }
protected:
    bcf_hdr_t *hdr() { return hdr_.get(); }

private:
    //std::shared_ptr<bcf_hdr_t, void(*)(bcf_hdr_t *)> hdr_;
    std::shared_ptr<bcf_hdr_t> hdr_;

public:
    // Indicates type of mutation in VCF record
    static const uint16_t REF = VCF_REF;
    static const uint16_t SNP = VCF_SNP;
    static const uint16_t MNP = VCF_MNP;
    static const uint16_t INDEL = VCF_INDEL;
    static const uint16_t OTHER = VCF_OTHER;

    friend class Variant;
};

inline Variant::Variant(const File &file) : BareVariant(), hdr_{file.hdr_} {
    bcf_float_set_missing(qual);
}


struct allele_t {
    operator int() {
        return value;
    }
    int value;
};

// Wrappers for htslib "genotype" functions
inline allele_t encode_allele_phased(int index) {
    return {bcf_gt_phased(index)};
}

inline allele_t encode_allele_unphased(int index) {
    return {bcf_gt_unphased(index)};
}

inline constexpr allele_t encode_allele_missing() {
    return {bcf_gt_missing};
}

inline bool allele_is_missing(allele_t value) {
    return bcf_gt_is_missing(value);
}

inline bool allele_is_phased(allele_t value) {
    return bcf_gt_is_phased(value);
}

inline int decode_allele(allele_t value) {
    return bcf_gt_allele(value);
}

inline int genotype_from_alleles(int a, int b) {
    return bcf_alleles2gt(a, b);
}

inline std::pair<int, int> alleles_from_genotype(int value) {
    int a, b;
    bcf_gt2alleles(value, &a, &b);
    return {a, b};
}

} // namespace bcf
} // namespace hts

#endif /* CXX_HTS_BCF_H */
