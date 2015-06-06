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

namespace hts {
namespace bcf {

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
    File(hts::File &&other) : hts::File(std::move(other)),
        hdr_{nullptr, bcf_hdr_destroy}, rec_{nullptr, bcf_destroy} {

        if(!is_open()) { // nothing to do
            return;
        }

        if(is_write()) {
            hdr_.reset(bcf_hdr_init("w"));
        } else {
            if(format().category != variant_data)
                throw std::runtime_error("file '" + std::string(name())
                                         + "' does not contain variant data (BCF/VCF).");
            hdr_.reset(bcf_hdr_read(handle()));
        }
        if(!hdr_) {
            throw std::runtime_error("unable to access header in file '" + std::string(
                                         name()) + "'.");
        }
        rec_.reset(bcf_init1());
        if(!rec_) {
            throw std::runtime_error("unable to construct bcf record.");
        }
    }

    File(const char *file, const char *mode) : File(hts::File(file, mode)) {
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
        bcf_hdr_add_sample(hdr(),
                           nullptr); // htslib requires NULL sample before it will write out all the other samples
        return bcf_hdr_write(handle(), hdr());
    }

    std::pair<char **, int> samples() const {
        return {hdr_->samples, bcf_hdr_nsamples(header())};
    }

    const bcf_hdr_t *header() const { return hdr_.get(); }

    //TODO: Split off rec into its own wrapper.

    /**
     * SetID() - Sets the CHROM, POS, and ID fields in the VCF record.
     * chrom is required field. If id is NULL then ID value will default to "."
     */
    void SetID(const char *chrom, int pos, const char *id = nullptr) {
        // CHROM
        if(chrom == nullptr) {
            return;    // should display some error?
        }
        rec_->rid = bcf_hdr_name2id(hdr(), chrom);

        // POS
        rec_->pos = pos;

        // ID
        if(id != nullptr) {
            bcf_update_id(hdr(), rec_.get(), id);
        }
    }

    void SetQuality(int quality) {
        rec_->qual = quality;
    }


    int SetFilter(const char *filter) {
        if(filter == nullptr) {
            return 0;
        }
        int32_t fid = bcf_hdr_id2int(hdr(), BCF_DT_ID, filter);
        return bcf_update_filter(hdr(), rec_.get(), &fid, 1);
    }

    /**
     * SetAlleles() - Set the REF and ALT fields in the current record
     * @str: A comma separated list of all the alleles that show up in the sample.
     *       the first element is the REF value.
     */
    int SetAlleles(const std::string &str) {
        if(str.empty()) {
            return 0;
        }
        return bcf_update_alleles_str(hdr(), rec_.get(), str.c_str());
    }


    /**
     * UpdateInfoField() - Add another key=Value pair to the INFO field
     * @key: must be defined in the Header or won't be added.
     * @value: if set to NULL then previous set key-value pairs will be removed.
     *
     * TODO: Add a vector<> version of flat, ints, and strings.
     */
    int UpdateInfo(const char *key, float value) {
        return bcf_update_info_float(hdr(), rec_.get(), key, &value, 1);
    }

    int UpdateInfo(const char *key, int32_t value) {
        return bcf_update_info_int32(hdr(), rec_.get(), key, &value, 1);
    }

    int UpdateInfo(const char *key, std::string &value) {
        return bcf_update_info_string(hdr(), rec_.get(), key, value.c_str());
    }


    /**
     * UpdateSample() - Add a key to the FORMAT field and update each sample column
     *                  with the corresponding key values.
     * @name: the tag name that appears in FORMAT field. Must be defined in the header
     * @data: A list of values that will populate the sample/genotype fields. The vector
     *        should be a multiple of the number sample columns.
     */
    int UpdateSamples(const char *name, const std::vector<float> &data) {
        assert(name != nullptr);
        return bcf_update_format_float(hdr(), rec_.get(), name, &data[0], data.size());
    }
    int UpdateSamples(const char *name, const std::vector<int32_t> &data) {
        assert(name != nullptr);
        return bcf_update_format_int32(hdr(), rec_.get(), name, &data[0], data.size());
    }
    int UpdateSamples(const char *name, const std::vector<const char *> &data) {
        assert(name != nullptr);
        return bcf_update_format_string(hdr(), rec_.get(), name,
                                        const_cast<const char **>(&data[0]), data.size());
    }
    int UpdateSamples(const char *name, const std::vector<std::string> &data) {
        std::vector<const char *> v(data.size());
        for(decltype(data.size()) u = 0; u < data.size(); ++u) {
            v[u] = data[u].c_str();
        }
        return UpdateSamples(name, v);
    }

    /** Writes out the up-to-date info in the record and prepares for the next line */
    void WriteRecord() {
        // Add line to the body of the VCF
        bcf_write1(handle(), hdr(), rec_.get());
        // reset the record for the next line
        bcf_clear(rec_.get());
    }
protected:
    bcf_hdr_t *hdr() { return hdr_.get(); }

private:
    std::unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t *)> hdr_;
    std::unique_ptr<bcf1_t, void(*)(bcf1_t *)> rec_;

public:
    // Indicates type of mutation in VCF record
    static const uint16_t REF = VCF_REF;
    static const uint16_t SNP = VCF_SNP;
    static const uint16_t MNP = VCF_MNP;
    static const uint16_t INDEL = VCF_INDEL;
    static const uint16_t OTHER = VCF_OTHER;

};

}
} // namespace hts::bcf

#endif /* CXX_HTS_BCF_H */
