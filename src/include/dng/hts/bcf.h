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
#include <cstdlib>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

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
const float float_vector_end = []() -> float {
    // We use a lambda function so the result can be constant.
    float a;
    bcf_float_set_vector_end(a);
    return a;
}();

constexpr int32_t int32_missing = bcf_int32_missing;
constexpr int32_t int32_vector_end = bcf_int32_vector_end;

constexpr int16_t int16_missing = bcf_int16_missing;
constexpr int16_t int16_vector_end = bcf_int16_vector_end;

constexpr int8_t int8_missing = bcf_int8_missing;
constexpr int8_t int8_vector_end = bcf_int8_vector_end;

const std::string str_missing = ".";

template<typename T>
using buffer_t = std::unique_ptr<T[],decltype(std::free)*>;

// bcf_get_format_* uses realloc internally, so this buffer
// will be managed by malloc and free
template<typename T>
buffer_t<T> make_buffer(std::size_t sz) {
    void *p = std::malloc(sizeof(T)*sz);
    if(p == nullptr) {
        throw std::bad_alloc{};
    }
    return { reinterpret_cast<T*>(p),std::free};
}

typedef bcf1_t BareVariant;

inline
int get_format_int32(const bcf_hdr_t *header, BareVariant *record, const char *tag, buffer_t<int>* buffer, int *capacity) {
    int *p = buffer->get();
    int n = bcf_get_format_int32(header, record, tag, &p, capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->get()) {
        // update pointer
        buffer->release();
        buffer->reset(p);
    }
    return n;
}

inline
int get_format_float(const bcf_hdr_t *header, BareVariant *record, const char *tag, buffer_t<float>* buffer, int *capacity) {
    float *p = buffer->get();
    int n = bcf_get_format_float(header, record, tag, &p, capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->get()) {
        // update pointer
        buffer->release();
        buffer->reset(p);
    }
    return n;
}

inline
int get_info_int32(const bcf_hdr_t *header, BareVariant *record, const char *tag, buffer_t<int>* buffer, int *capacity) {
    int *p = buffer->get();
    int n = bcf_get_info_int32(header, record, tag, &p, capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->get()) {
        // update pointer
        buffer->release();
        buffer->reset(p);
    }
    return n;
}

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

    // Turn these off until we figure out how to use bcf_copy
    BareVariant& operator=(const BareVariant&) = delete;
    Variant(const BareVariant&) = delete;

    // Getters
    int32_t target_id() const { return rid; }
    int32_t position() const { return pos; }
    int32_t ref_length() const { return rlen; }
    float quality() const { return qual; }

    int num_alleles() const {
        return base()->n_allele;
    }
    const char* allele(int n) const {
        assert(n < num_alleles());
        
        return base()->d.allele[n];
    }

    // Setters
    void quality(float q) { qual = q; }
    void target(const char *chrom) {
        rid = (chrom != nullptr) ? bcf_hdr_name2id(header(), chrom) : -1;
    }
    void position(int32_t p) { pos = p; }
    void id(const char *str) {
        bcf_update_id(header(), base(), str);
    }
    int filter(const char *str) {
        if(str == nullptr) {
            return 0;
        }
        int32_t fid = bcf_hdr_id2int(header(), BCF_DT_ID, str);
        return bcf_add_filter(header(), base(), fid);
    }
    int filter(const std::string &str) {
        int32_t fid = bcf_hdr_id2int(header(), BCF_DT_ID, str.c_str());
        return bcf_add_filter(header(), base(), fid);
    }

    /**
     * alleles() - Set the REF and ALT fields in the current record
     * @str: A comma separated list of all the alleles that show up in the sample.
     *       the first element is the REF value.
     */
    int alleles(const std::string &str) {
        return bcf_update_alleles_str(header(), base(), str.c_str());
    }
    int alleles(const char *str) {
        if(str == nullptr) {
            return 0;
        }
        return bcf_update_alleles_str(header(), base(), str);
    }
    int alleles(const char **alleles, int num_alleles) {
        if(alleles == nullptr) {
            return 0;
        }
        return bcf_update_alleles(header(), base(), alleles, num_alleles);
    }

    /**
     * info() - Add another key=Value pair to the INFO field
     * @key: must be defined in the Header or won't be added.
     * @value: if set to NULL then previous set key-value pairs will be removed.
     */
    int info(const char *key, const float *value, std::size_t count) {
        assert(key != nullptr);
        return bcf_update_info_float(header(), base(), key, value, count);
    }
    int info(const char *key, const int32_t *value, std::size_t count) {
        assert(key != nullptr);
        return bcf_update_info_int32(header(), base(), key, value, count);
    }
    int info(const char *key, const std::string &value) {
        assert(key != nullptr);
        return bcf_update_info_string(header(), base(), key, value.c_str());
    }
    int info(const char *key, const char *value) {
        assert(key != nullptr);
        return bcf_update_info_string(header(), base(), key, value);
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
    int samples(const char *name, const float* data, size_t sz) {
        assert(name != nullptr);
        return bcf_update_format_float(header(), base(), name, data, sz);
    }
    int samples(const char *name, const int32_t* data, size_t sz) {
        assert(name != nullptr);
        return bcf_update_format_int32(header(), base(), name, data, sz);
    }
    int samples(const char *name, const char* const* data, size_t sz) {
        assert(name != nullptr);
        return bcf_update_format_string(header(), base(), name,
            const_cast<const char **>(data), sz);
    }    
    template<typename T, typename A>
    int samples(const char *name, const std::vector<T, A> &data) {
        return samples(name, data.data(), data.size());
    }
    template<typename T, size_t N>
    int samples(const char *name, const std::array<T, N> &data) {
        return samples(name, data.data(), data.size());
    }
    template<typename T, size_t N>
    int samples(const char *name, const T(&data)[N]) {
        return samples(name, &data[0], N);
    }    

    int samples(const char *name, const std::vector<std::string> &data) {
        std::vector<const char *> v(data.size());
        for(decltype(data.size()) u = 0; u < data.size(); ++u) {
            v[u] = data[u].c_str();
        }
        return samples(name, v);
    }
    int sample_genotypes(const std::vector<int32_t> &data) {
        return bcf_update_genotypes(header(), base(), &data[0], data.size());
    }

protected:
    BareVariant *base() {return static_cast<BareVariant *>(this);}
    const BareVariant *base() const {return static_cast<const BareVariant *>(this);}

    const bcf_hdr_t *header() {
        // check if the header pointer has not been initialized
        assert(header_);
        return header_.get();
    }

private:
    std::shared_ptr<const bcf_hdr_t> header_;

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
            header_ = std::shared_ptr<bcf_hdr_t>(bcf_hdr_init("w"), bcf_hdr_destroy);
        } else {
            if(format().category != variant_data)
                throw std::runtime_error("file '" + std::string(name())
                                         + "' does not contain variant data (BCF/VCF).");
            header_ = std::shared_ptr<bcf_hdr_t>(bcf_hdr_read(handle()), bcf_hdr_destroy);
        }
        if(!header_) {
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
        return bcf_hdr_append(header(), line);
    }
    int AddHeaderMetadata(const std::string &line) {
        return bcf_hdr_append(header(), line.c_str());
    }

    /** AddHeaderMetadata() - Adds a "##key=value" line to the VCF header */
    int AddHeaderMetadata(const char *key, const char *value) {
        if(key == nullptr || value == nullptr) {
            return -1;
        }
        std::string line = std::string("##") + key + "=" + value;
        return bcf_hdr_append(header(), line.c_str());
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
        return bcf_hdr_add_sample(header(), sample);
    }

    /** Add a "#contig=" metadata entry to the VCF header. */
    int AddContig(const char *contigid, uint32_t length) {
        if(contigid == nullptr) {
            return -1;
        }
        std::string conv = std::string("##contig=<ID=") + contigid
                           + ",length=" + std::to_string(length) + ">";
        return bcf_hdr_append(header(), conv.c_str());
    }

    /** Creates the header field. Call only after adding all the sample fields */
    int WriteHeader() {
        // htslib requires NULL sample before it will write out all the other samples
        bcf_hdr_add_sample(header(), nullptr);
        return bcf_hdr_write(handle(), header());
    }

    std::pair<char **, int> samples() const {
        return {header()->samples, bcf_hdr_nsamples(header())};
    }

    std::vector<std::pair<const char *, int>> contigs() const;

    const bcf_hdr_t *header() const { return header_.get(); }

    /** Writes out the up-to-date info in the record and prepares for the next line */
    void WriteRecord(Variant &rec) {
        // Add line to the body of the VCF
        assert(rec.header() == header());
        bcf_write(handle(), header(), &rec);
    }

    void WriteRecord(BareVariant *rec) {
        // Add line to the body of the VCF
        assert(rec != nullptr);
        bcf_write(handle(), header(), rec);
    }

    /** Explicity closes the file and flushes the stream */
    void Close() {
    	bcf_close(handle());
    }


protected:
    bcf_hdr_t *header() { return header_.get(); }

private:
    //std::shared_ptr<bcf_hdr_t, void(*)(bcf_hdr_t *)> hdr_;
    std::shared_ptr<bcf_hdr_t> header_;

public:
    // Indicates type of mutation in VCF record
    static const uint16_t REF = VCF_REF;
    static const uint16_t SNP = VCF_SNP;
    static const uint16_t MNP = VCF_MNP;
    static const uint16_t INDEL = VCF_INDEL;
    static const uint16_t OTHER = VCF_OTHER;

    friend class Variant;
};

inline
std::vector<std::pair<const char *, int>> contigs(const bcf_hdr_t *header) {
    assert(header != nullptr);
    const int num_contigs = header->n[BCF_DT_CTG];

    std::vector<std::pair<const char *, int>> contigs;

    for(int i=0;i<num_contigs;++i) {
        const char *name = header->id[BCF_DT_CTG][i].key;
        int length = header->id[BCF_DT_CTG][i].val->info[0];
        contigs.emplace_back(name, length);
    }
    return contigs;
}

inline
std::vector<std::pair<const char *, int>> File::contigs() const {
    return bcf::contigs(header());
}


inline
Variant::Variant(const File &file) : BareVariant(), header_{file.header_} {
    bcf_float_set_missing(qual);
}

class SyncedReader {
public:
    enum struct Collapse {
        None = COLLAPSE_NONE,
        SNPs = COLLAPSE_SNPS,
        Indels = COLLAPSE_INDELS,
        Any = COLLAPSE_ANY,
        Some = COLLAPSE_SOME,
        Both = COLLAPSE_BOTH
    };

    SyncedReader(Collapse collapse = Collapse::None) : handle_{bcf_sr_init(), bcf_sr_destroy} {
        handle()->collapse = static_cast<int>(collapse);
    }

    bcf_srs_t * handle() {
        assert(handle_);
        return handle_.get();
    }
    const bcf_srs_t * handle() const {
        assert(handle_);
        return handle_.get();
    }

    int SetRegions(const char* regions) {
        assert(regions != nullptr);
        assert(num_readers() == 0); // must be called before any readers are added
        return bcf_sr_set_regions(handle(), regions, 0);
    }

    int AddReader(const char *filename) {
        assert(filename != nullptr);
        return bcf_sr_add_reader(handle(), filename);
    }
    void RemoveReader(int index) {
        assert(0 <= index && index < handle()->nreaders);
        bcf_sr_remove_reader(handle(), index);
    }

    int NextLine() {
        return bcf_sr_next_line(handle());
    }
    bcf1_t* GetLine(int index) {
        assert(0 <= index && index < handle()->nreaders);
        return bcf_sr_get_line(handle(), index);
    }

    int SetSamples(const char* samples) {
        assert(samples != nullptr);
        return bcf_sr_set_samples(handle(), samples, 0);
    }

    int num_readers() const {
        return handle()->nreaders;
    }    
    bcf_sr_t* reader(int index) {
        assert(0 <= index && index < handle()->nreaders);
        return bcf_sr_get_reader(handle(), index);    
    }
    const bcf_sr_t* reader(int index) const {
        assert(0 <= index && index < handle()->nreaders);
        return bcf_sr_get_reader(handle(), index);    
    }
    bcf_hdr_t* header(int index) {
        assert(0 <= index && index < handle()->nreaders);
        return bcf_sr_get_header(handle(), index);    
    }
    const bcf_hdr_t* header(int index) const {
        assert(0 <= index && index < handle()->nreaders);
        return bcf_sr_get_header(handle(), index);    
    }

    const char * const * samples() const {
        return handle()->samples;
    }
    int num_samples() const {
        return handle()->n_smpl;
    }


private:
    std::unique_ptr<bcf_srs_t, void(*)(bcf_srs_t *)> handle_;
};

// Use a structure to have strongly-typed alleles
struct allele_t {
    constexpr operator int() const {
        return value;
    }
    int value;
};

// Wrappers for htslib "genotype" functions
constexpr inline allele_t encode_allele_phased(int index) {
    return {bcf_gt_phased(index)};
}

constexpr inline allele_t encode_allele_unphased(int index) {
    return {bcf_gt_unphased(index)};
}

constexpr inline allele_t encode_allele_missing() {
    return {bcf_gt_missing};
}

constexpr inline bool allele_is_missing(allele_t value) {
    return bcf_gt_is_missing(value);
}

constexpr inline bool allele_is_phased(allele_t value) {
    return bcf_gt_is_phased(value);
}

constexpr inline int decode_allele(allele_t value) {
    return bcf_gt_allele(value);
}

constexpr inline int genotype_from_alleles(int a, int b) {
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
