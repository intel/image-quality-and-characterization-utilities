/*
* Copyright (c) 2019, Intel Corporation
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the Intel Corporation nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#pragma once

#include "Teisko/Image/API.hpp"
#include <cstdint>      // uint16_t, uint32_t etc
#include <cstring>      // memcpy
#include <iterator>
#include <string>       // std::string
#include <memory>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <exception>

namespace Teisko
{
    /// Some encountered tiff tags
    enum tiff_tags
    {
        tag_new_subfile_type = 0x00fe,
        tag_image_width = 0x0100,
        tag_image_height = 0x0101,
        tag_bits_per_sample = 0x102,
        tag_compression_method = 0x103,
        /**/    tag_compression_method_uncompressed = 1,
        tag_photometric_interpretation = 0x106,
        /**/    tag_photometric_interpretation_bayer = 1,
        /**/    tag_photometric_interpretation_rgb = 2,
        /**/    tag_photometric_interpretation_yuv = 6,
        tag_image_description = 0x10e,              // ASCII - support not ready
        tag_strip_offsets = 0x0111,
        tag_orientation = 0x0112,
        /**/    tag_orientation_top_left = 1,
        tag_samples_per_pixel = 0x0115,
        tag_rows_per_strip = 0x0116,                // set to image height (or -1 == infinite)
        tag_strip_byte_counts = 0x0117,             // image size in bytes
        tag_minimum_sample_value = 0x0118,
        tag_maximum_sample_value = 0x0119,
        tag_x_resolution = 0x011a,
        tag_y_resolution = 0x011b,
        tag_planar_configuration = 0x011c,
        /**/    tag_planar_configuration_chunky = 1,    // single plane
        /**/    tag_planar_configuration_planar = 2,    // multiple planes (not supported)
        tag_resolution_unit = 0x0128,
        tag_software = 0x131,
        tag_date_and_time = 0x132,
        tag_predictor = 0x013d,
        tag_tile_offsets = 0x0144,
        tag_tile_byte_counts = 0x0145,
        tag_sample_format = 0x153,
        /**/    tag_sample_format_unsigned = 1,
        /// Minimal set of extended tags required for YCbCr file
        tag_ycbcr_coefficients = 0x211,           // rational: 299 / 1000, 587 / 1000, 114 / 1000
        tag_ycbcr_subsampling = 0x212,            // array of shorts : [1 1]
        tag_ycbcr_positioning = 0x213,
        /**/    tag_ycbcr_positioning_center = 1,
        tag_reference_blackwhite = 0x214,         // 6 rationals
        tag_xmp = 0x02bc,
        tag_iptc_metadata = 0x83bb,             // Undefined or Byte
        tag_photoshop = 0x8649,
        tag_exif_ifd = 0x8769,                  // Offset to Exif IFD: [count][tags * M][next]
        tag_icc_profile = 0x8773,
        tag_exif_colorspace = 0xa001,
        tag_exif_pixel_x_dimension = 0xa002,
        tag_exif_pixel_y_dimension = 0xa003
    };

    enum tiff_tag_types
    {
        BYTE = 1,
        ASCII = 2,
        SHORT,
        LONG,
        RATIONAL,
        SBYTE,
        UNDEFINED,   // Binary with custom meaning
        SSHORT,
        SLONG,
        SRATIONAL,
        FLOAT,
        DOUBLE
    };

    /// Template class to associate some real object to each enumerated type
    template<tiff_tag_types> struct tiff_type_;

    // specialization for each type
    template<> struct tiff_type_<BYTE> { using type = uint8_t; };
    template<> struct tiff_type_<ASCII> { using type = char; };
    template<> struct tiff_type_<SHORT> { using type = uint16_t; };
    template<> struct tiff_type_<LONG> { using type = uint32_t; };
    template<> struct tiff_type_<RATIONAL> { using type = uint32_t; };
    template<> struct tiff_type_<SBYTE> { using type = signed char; };
    template<> struct tiff_type_<UNDEFINED> { using type = uint8_t; };
    template<> struct tiff_type_<SSHORT> { using type = int16_t; };
    template<> struct tiff_type_<SLONG> { using type = int32_t; };
    template<> struct tiff_type_<SRATIONAL> { using type = int32_t; };
    template<> struct tiff_type_<FLOAT> { using type = float; };
    template<> struct tiff_type_<DOUBLE> { using type = double; };

    const int tiff_tag_size = 12;
    /// Check if given byte order matches native order
    /// Used to decide if input data has to be byte swapped
    static bool is_native(bool is_bigendian)
    {
        union {
            uint16_t w;
            uint8_t b[2];
        } temp;
        temp.w = 42;
        return is_bigendian == (temp.b[1] == 42);
    }

    /// Extracts 16-bit word to native format
    static uint16_t get_short(uint8_t *data, bool is_bigendian)
    {
        return is_bigendian ?
            ((uint16_t)data[0] << 8) | ((uint16_t)data[1] << 0) :
            ((uint16_t)data[0] << 0) | ((uint16_t)data[1] << 8);
    }

    /// Extracts 32-bit word to native format
    static uint32_t get_int(uint8_t *data, bool is_bigendian)
    {
        return is_bigendian ?
            ((uint32_t)data[0] << 24) | ((uint32_t)data[1] << 16) |
            ((uint32_t)data[2] << 8) | ((uint32_t)data[3] << 0) :
            ((uint32_t)data[3] << 24) | ((uint32_t)data[2] << 16) |
            ((uint32_t)data[1] << 8) | ((uint32_t)data[0] << 0);
    }
    /// Writes 16-bit word in network/little endian format
    static void put_short(std::ostream &stream, uint16_t data, bool is_bigendian)
    {
        if (is_bigendian)
            stream.put(data >> 8).put(data & 0xff);
        else
            stream.put(data & 0xff).put(data >> 8);
    }
    /// Writes 32-bit word in network/little endian format
    static void put_int(std::ostream &out, uint32_t data, bool is_bigendian)
    {
        if (is_bigendian)
            out.put(data >> 24).put((data >> 16) & 0xff).put((data >> 8) & 0xff).put(data & 0xff);
        else
            out.put(data & 0xff).put((data >> 8) & 0xff).put((data >> 16) & 0xff).put(data >> 24);
    }

    /// container class for all tags in Tiff file
    /// Contains ID, type and value
    struct tiff_tag_base
    {
        /// Default constructor
        tiff_tag_base(
            tiff_tags id = (tiff_tags)0,
            tiff_tag_types type = (tiff_tag_types)0,
            uint32_t count = 0,
            uint32_t unit_size = 0) :
            _id((uint16_t)id),
            _type((uint16_t)type),
            _count(count),
            _offset(0),
            _swap_length(unit_size)
        { }

        /// Construct a tag from input stream
        tiff_tag_base(std::vector<uint8_t> stream, uint32_t offset, bool is_bigendian = false)
        {
            if (offset + tiff_tag_size > stream.size())
                throw std::runtime_error("Offset points outside of tiff data");

            uint8_t *dptr = stream.data();

            _id = get_short(dptr + offset, is_bigendian);
            _type = get_short(dptr + offset + 2, is_bigendian);
            _count = get_int(dptr + offset + 4, is_bigendian);
            _offset = 0;

            // Validate type and get length of elementary unit
            _swap_length = element_size(_type);

            uint32_t total_size = _swap_length * _count;
            total_size <<= (_type == RATIONAL || _type == SRATIONAL) ? 1 : 0;

            offset += 8;
            if (total_size > 4)
                offset = get_int(dptr + offset, is_bigendian);

            if (offset + total_size > stream.size())
                throw std::runtime_error("Offset to values points outside of tiff data");

            uint32_t swap_mask = is_native(is_bigendian) ? 0 : _swap_length - 1;
            _additional_data = std::vector<uint8_t>(total_size);
            uint8_t *dptr_dst = _additional_data.data();

            // copy the string byte per byte, swapping byte order when necessary
            for (uint32_t i = 0; i < total_size; i++)
                dptr_dst[i] = dptr[offset + (i ^ swap_mask)];
        }

        /// Return tag value as string
        inline std::string to_string();

        /// Return tag value as vector of chosen type
        template <typename T = int> std::vector<T> value();

        /// Return tag ID
        tiff_tags tag() const { return (tiff_tags)_id; }

        /// Places tag data into tiff section pointed by 'offset'
        /// if the data doesn't fit in the fixed length tag structure
        void fix_offset(uint32_t &offset)
        {
            uint32_t sz = (uint32_t)_additional_data.size();
            if (sz > 4)
            {
                _offset = offset;
                offset += ((sz + 1) >> 1) << 1;
            }
        }

        /// Writes data that doesn't fit in the fixed length tag structure
        void write_additional_data(std::ostream &stream, bool is_bigendian)
        {
            uint32_t sz = (uint32_t)_additional_data.size();
            if (sz <= 4)
                return;

            uint8_t *dptr = _additional_data.data();
            uint32_t swap_mask = is_native(is_bigendian) ? 0 : _swap_length - 1;

            for (uint32_t i = 0; i < sz; i++)
                stream.put(dptr[i ^ swap_mask]);

            if (sz & 1)
                stream.put(0);      // Pad to word boundary
        }

        // Writes ID/TYPE/COUNT and OFFSET to output stream
        //  - populates OFFSET by data, if it fits
        void write_pod(std::ostream &stream, bool is_network_order)
        {
            put_short(stream, _id, is_network_order);
            put_short(stream, _type, is_network_order);
            put_int(stream, _count, is_network_order);
            uint32_t sz = (uint32_t)_additional_data.size();
            if (sz <= 4)
            {
                // Write up to 4 bytes from the data vector; then pad to 4-byte boundary
                uint32_t swap_mask = is_native(is_network_order) ? 0 : _swap_length - 1;
                for (uint32_t i = 0; i < 4; i++)
                    stream.put(i < sz ? _additional_data[i ^ swap_mask] : 0);
            }
            else
                put_int(stream, _offset, is_network_order);
        }

    protected:
        uint16_t _id;
        uint16_t _type;
        uint32_t _count;
        uint32_t _offset;
        uint32_t _swap_length;
        std::vector<uint8_t> _additional_data;

        // Validates type -- throws on error
        // Returns size of elementary unit
        inline static int element_size(uint16_t type);
    };

    // We can't achieve real run time polymorphism with this approach
    // but specializing the tag allows parameter type checking / error highlighting
    // - e.g. tiff_tag<SHORT>(tag_image_width, { 123456 });  is checked by IDE
    //    - the value 123456 doesn't fit the (unsigned short) field
    template <tiff_tag_types tag_type>
    struct tiff_tag : public tiff_tag_base
    {
        // Disallow (at compile time) specializations outside range [1..12]
        static_assert((int)tag_type > 0 && (int)tag_type <= (int)DOUBLE,
            "Tag types between 1 and 12 supported");

        using enum_type = typename tiff_type_<tag_type>::type;

        // Populate the values as iterator to [first, last[
        // All other variants are converted to this version
        // RATIONAL types strip odd number
        tiff_tag(tiff_tags id, const enum_type *first, const enum_type *last) :
            tiff_tag_base(id, tag_type, (uint32_t)(last - first), (uint32_t)sizeof(enum_type))
        {
            uint32_t shift = (tag_type == RATIONAL || tag_type == SRATIONAL) ? 1 : 0;
            _count >>= shift;       // Truncating divide by two (or nothing)
            _additional_data = std::vector<uint8_t>(
                reinterpret_cast<const uint8_t*>(first),
                reinterpret_cast<const uint8_t*>(first + (_count << shift)));
        }

        // Populate the values as iterator to first, length of array
        tiff_tag(tiff_tags id, const enum_type *first, size_t length) :
            tiff_tag(id, first, first + length) { }

        // Populate the values as initializer list -- useful for few values
        tiff_tag(tiff_tags id, std::initializer_list<enum_type> list) :
            tiff_tag(id, list.begin(), list.end()) { }

        // Populate tag from vector
        tiff_tag(tiff_tags id, const std::vector<enum_type> &vec) :
            tiff_tag(id, vec.data(), vec.data() + vec.size()) { }

        // Populate tag by single value
        tiff_tag(tiff_tags id, const enum_type &value) :
            tiff_tag(id, &value, &value + 1) { }

        // Populate tag by string -- not applicable to all types?
        tiff_tag(tiff_tags id, const std::string &str) :
            tiff_tag(id, str.c_str(), str.c_str() + str.length() + 1) { }

        // Converts all array / single value types to strings
        // ASCII has a specialization
        std::string to_string()
        {
            std::string result = "";
            size_t size = _additional_data.size() / _swap_length;
            auto begin = reinterpret_cast<const enum_type*>(_additional_data.data());
            char delimiter = (_type == RATIONAL || _type == SRATIONAL) ? '/' : ' ';
            char toggle = delimiter ^ ' ';  // toggles between "/" and " " or ' ' and ' '
            for (size_t i = 0; i < size; i++)
            {
                result += std::to_string(*begin++);
                if (i != size - 1)
                    result += delimiter;
                delimiter ^= toggle;
            }
            return result;
        }

        // Converts the data to any type T allowing e.g. narrowing conversion
        template <typename T = enum_type>
        std::vector<T> value()
        {
            size_t size = _additional_data.size() / sizeof(enum_type);
            auto begin = reinterpret_cast<const enum_type*>(_additional_data.data());
            // Iterate the vector element by element to allow narrowing conversion
            // alternative:  result std::vector<T>(begin, end)  would not compile for all types
            std::vector<T> result(size);
            std::generate(result.begin(), result.end(), [&begin]() { return static_cast<T>(*begin++); });
            return result;
        }
    };

    /// Specialization for ASCII -- convert buffer to string (without trailing ASCIIZ)
    template <>
    inline std::string tiff_tag<ASCII>::to_string() {
        size_t len = _additional_data.size();
        auto begin = _additional_data.begin();
        return len > 0 ? std::string(begin, begin + len - 1) : "";
    };

    template <typename T>
    std::vector<T> tiff_tag_base::value()
    {
        switch (_type)
        {
        case BYTE: return reinterpret_cast<tiff_tag<BYTE>*>(this)->value<T>();
        case ASCII: return reinterpret_cast<tiff_tag<ASCII>*>(this)->value<T>();
        case SHORT: return reinterpret_cast<tiff_tag<SHORT>*>(this)->value<T>();
        case LONG: return reinterpret_cast<tiff_tag<LONG>*>(this)->value<T>();
        case RATIONAL: return reinterpret_cast<tiff_tag<RATIONAL>*>(this)->value<T>();
        case SBYTE: return reinterpret_cast<tiff_tag<SBYTE>*>(this)->value<T>();
        case UNDEFINED: return reinterpret_cast<tiff_tag<UNDEFINED>*>(this)->value<T>();
        case SSHORT: return reinterpret_cast<tiff_tag<SSHORT>*>(this)->value<T>();
        case SLONG: return reinterpret_cast<tiff_tag<SLONG>*>(this)->value<T>();
        case SRATIONAL: return reinterpret_cast<tiff_tag<SRATIONAL>*>(this)->value<T>();
        case FLOAT: return reinterpret_cast<tiff_tag<FLOAT>*>(this)->value<T>();
        case DOUBLE: return reinterpret_cast<tiff_tag<DOUBLE>*>(this)->value<T>();
        default:
            throw std::runtime_error("Invalid tiff type");
        }
    };

    std::string tiff_tag_base::to_string()
    {
        switch (_type)
        {
        case BYTE: return reinterpret_cast<tiff_tag<BYTE>*>(this)->to_string();
        case ASCII: return reinterpret_cast<tiff_tag<ASCII>*>(this)->to_string();
        case SHORT: return reinterpret_cast<tiff_tag<SHORT>*>(this)->to_string();
        case LONG: return reinterpret_cast<tiff_tag<LONG>*>(this)->to_string();
        case RATIONAL: return reinterpret_cast<tiff_tag<RATIONAL>*>(this)->to_string();
        case SBYTE: return reinterpret_cast<tiff_tag<SBYTE>*>(this)->to_string();
        case UNDEFINED: return reinterpret_cast<tiff_tag<UNDEFINED>*>(this)->to_string();
        case SSHORT: return reinterpret_cast<tiff_tag<SSHORT>*>(this)->to_string();
        case SLONG: return reinterpret_cast<tiff_tag<SLONG>*>(this)->to_string();
        case SRATIONAL: return reinterpret_cast<tiff_tag<SRATIONAL>*>(this)->to_string();
        case FLOAT: return reinterpret_cast<tiff_tag<FLOAT>*>(this)->to_string();
        case DOUBLE: return reinterpret_cast<tiff_tag<DOUBLE>*>(this)->to_string();
        default:
            throw std::runtime_error("Invalid tiff type");
        }
    }

    int tiff_tag_base::element_size(uint16_t type)
    {
        switch (type)
        {
        case BYTE: return sizeof(tiff_type_<BYTE>::type);
        case ASCII: return sizeof(tiff_type_<ASCII>::type);
        case SHORT: return sizeof(tiff_type_<SHORT>::type);
        case LONG: return sizeof(tiff_type_<LONG>::type);
        case RATIONAL: return sizeof(tiff_type_<RATIONAL>::type);
        case SBYTE: return sizeof(tiff_type_<SBYTE>::type);
        case UNDEFINED: return sizeof(tiff_type_<UNDEFINED>::type);
        case SSHORT: return sizeof(tiff_type_<SSHORT>::type);
        case SLONG: return sizeof(tiff_type_<SLONG>::type);
        case SRATIONAL: return sizeof(tiff_type_<SRATIONAL>::type);
        case FLOAT: return sizeof(tiff_type_<FLOAT>::type);
        case DOUBLE: return sizeof(tiff_type_<DOUBLE>::type);
        default:
            throw std::runtime_error("Invalid tiff type");
        }
    }

    /// Container for all data found from tiff file
    /// The overall structure of tiff file (with single image) is
    ///         +-----+-----+-------------------+
    /// HEADER  | II  | 42  |offset to IFD #0   |
    ///         +-----+-----+-------------------+
    ///
    ///         +-----------------//---+
    /// Image   | IMAGE DATA           |
    ///         +-----------------//---+
    ///
    /// [EXIF   +---------------------------+  +---+---//------------+-----+
    ///  DATA]  | DATA pool for Exif Values |  | N | N * IFD_ENTRIES |  0  |
    ///         +---------------------------+  +---+-----------------+-----+
    ///
    /// Tags    +---------------//-------+
    ///         | Data pool for tag data |
    ///         +---------------//-------+
    ///
    /// IFD #0  +-------+-------//---------+------------------+
    ///         | Count | M * IFD_ENTRIES  | offset to IFD #1 |
    ///         +-------+-------//---------+------------------+
    ///
    ///  Specifically -- tag_exif_ifd.value == offset to EXIF_DATA
    ///               -- tag_strip_offsets == offset to Image data
    struct tiff_file
    {
        /// Make a tiff header for reading/writing
        explicit tiff_file(bool use_big_endian = false) : is_network_order(use_big_endian) { }

        /// Adds a tag to tiff_file
        void add_tag(tiff_tag_base &&tag)
        {
            _ifd[tag.tag()] = tag;
        }
        /// Adds multiple tags to tiff_file
        void add_tag(std::list<tiff_tag_base> &&tags)
        {
            for (auto &t : tags)
                add_tag(std::move(t));
        }
        /// Removes tag from directory
        void erase_tag(tiff_tags id)
        {
            auto it = _ifd.find(id);
            if (it != _ifd.end())
                _ifd.erase(it);
        }

        /// Returns a vector of all tags related to primary image
        std::vector<tiff_tags> tags()
        {
            std::vector<tiff_tags> tags;
            tags.reserve(_ifd.size());
            std::transform(_ifd.begin(), _ifd.end(), std::back_inserter(tags),
                [](const std::pair<tiff_tags, tiff_tag_base> &key_value) { return key_value.first; });
            return tags;
        }

        /// Returns the content of a tag as a string
        std::string tag_as_string(tiff_tags id)
        {
            auto it = _ifd.find(id);
            if (it == _ifd.end())
                return "";
            return it->second.to_string();
        }

        /// Returns tag values as a vector -- internal usage
        std::vector<int> tag_as_vector(tiff_tags id)
        {
            auto it = _ifd.find(id);
            if (it == _ifd.end())
                return std::vector<int>();
            else
                return it->second.value();
        }

        /// \brief  write_image   Writes an rgb or yuv image to tiff file
        /// \param  out           Byte serializer (e.g. opened file handle) for output
        /// \param  image         Array of image_api color planes
        /// \param  bpp           Original bits per sample
        /// \param  type          RGB or YUV
        /// The image will be pre-scaled from bpp to next bit depth in [8,16,32]
        /// Throws if all planes are not of equal size
        template <typename T>
        void write_image(
            std::ostream &out,
            image<T>(&input)[3],
            int bpp = sizeof(T) * 8,
            int type = tag_photometric_interpretation_rgb)
        {
            auto width = input[0]._width;
            auto height = input[0]._height;
            uint16_t hscale = 1;
            uint16_t vscale = 1;

            if (type == tag_photometric_interpretation_rgb)
            {
                if (width != input[1]._width ||
                    height != input[1]._height ||
                    width != input[2]._width ||
                    height != input[2]._height)
                    throw std::invalid_argument("All RGB image planes must match in dimensions");
            }
            else if (type == tag_photometric_interpretation_yuv)
            {
                if (input[1]._width != input[2]._width ||
                    input[1]._height != input[2]._height)
                {
                    throw std::invalid_argument("Planes 1 and 2 must match in dimensions");
                }
                hscale = check_factor(width, input[1]._width);
                vscale = check_factor(height, input[1]._height);

                add_tag(tiff_tag<SHORT>(tag_ycbcr_subsampling, { hscale, vscale }));

            }
            else
            {
                throw std::invalid_argument("Planar format only supported for rgb and yuv");
            }

            // Serialize the image data to 8,16 or 32 bits
            if (bpp <= 8)
            {
                serialize_image<uint8_t, T>(input, bpp, hscale, vscale);
                write_buffer(out, width, height, 8, type);
            }
            else if (bpp <= 16)
            {
                serialize_image<uint16_t, T>(input, bpp, hscale, vscale);
                write_buffer(out, width, height, 16, type);
            }
            else
            {
                serialize_image<uint32_t, T>(input, bpp, hscale, vscale);
                write_buffer(out, width, height, 32, type);
            }
        }

        /// \brief  write_image   Writes an rgb, yuv or Bayer image to tiff file
        /// \param  out_stream    Byte serializer (e.g. opened file handle) for output
        /// \param  image         Array of image_api color planes
        /// \param  bpp           Original bits per sample
        /// \param  type          RGB or YUV
        /// The image will be pre-scaled from bpp to next bit depth in [8,16,32]
        /// Throws if all planes are not of equal size
        template <typename T>
        void write_image(
            std::ostream &out,
            image<T> &input,
            int bpp = sizeof(T) * 8,
            int type = tag_photometric_interpretation_rgb)
        {
            int width = input._width;
            int height = input._height;

            if (type == tag_photometric_interpretation_rgb ||
                type == tag_photometric_interpretation_yuv)
            {
                if (input._width % 3 != 0)
                    throw std::invalid_argument("Scan line width must be multiple of three");
                width /= 3;
            }
            else if (type != tag_photometric_interpretation_bayer)
            {
                throw std::invalid_argument("Only supporting RGB, Bayer and YCbCr formats");
            }

            // Serialize the image data to 8,16 or 32 bits
            if (bpp <= 8)
            {
                serialize_image<uint8_t, T>(input, bpp);
                write_buffer(out, width, height, 8, type);
            }
            else if (bpp <= 16)
            {
                serialize_image<uint16_t, T>(input, bpp);
                write_buffer(out, width, height, 16, type);
            }
            else
            {
                serialize_image<uint32_t, T>(input, bpp);
                write_buffer(out, width, height, 32, type);
            }
        }

        /// Deserializes tiff stream to internal data
        /// Supports reading of first IFD, first exif IFD and
        /// image data using any compression method; tiles or strips
        void read(std::istream &input)
        {
            _ifd.clear();
            _exif.clear();
            _strips.clear();

            auto data = read_vector(input);

            auto offset = read_ifd_offset(data);
            if (offset > 0)
            {
                _ifd = read_ifd(data, offset);

                auto exif = tag_as_vector(tag_exif_ifd);

                if (exif.size() == 1)
                {
                    uint32_t offset_to_exif_ifd = exif[0];
                    _exif = read_ifd(data, offset_to_exif_ifd);
                }

                // First try to read image as strips
                if (read_image_data(data, tag_strip_byte_counts, tag_strip_offsets) == false)
                {
                    // And the try to read image as tiles
                    read_image_data(data, tag_tile_byte_counts, tag_tile_offsets);
                }
            }
        }

        /// Writes tiff tags and data to a stream
        void write(std::ostream &out_stream)
        {
            uint32_t offset = 8;
            tiff_make_image_data(offset);

            rebuild_exif_data(offset);

            fix_tag_offsets(_ifd, offset);      // No new tags should come after this...

            write_tiff_header(out_stream, offset);

            write_image_strips(out_stream);

            if (_exif.size() > 0)
            {
                write_additional_data(out_stream, _exif);
                write_ifd(out_stream, _exif);
            }

            write_additional_data(out_stream, _ifd);

            write_ifd(out_stream, _ifd);
        }

    private:
        bool is_network_order = false;

        std::map<tiff_tags, tiff_tag_base> _ifd;          // Only one image supported
        std::map<tiff_tags, tiff_tag_base> _exif;         // Only on exif structure supported
        std::vector<std::vector<uint8_t>> _strips;   // Serialized image -- any compression scheme

        /// \brief read_vector      Reads everything from an input stream
        /// \param  input           Input stream (e.g. memorystream, opened file)
        /// \returns                Vector of uint8_t
        static std::vector<uint8_t> read_vector(std::istream &input)
        {
            input.unsetf(std::ios::skipws);

            input.seekg(0, std::ios::end);
            std::streampos file_len = input.tellg();
            input.seekg(0, std::ios::beg);
            // Don't handle at all files larger than 4GB (no big-tiff support)
            if (file_len < 8 || file_len >= 0x100000000L)
                file_len = 0;
            std::vector<uint8_t> data((uint32_t)file_len);
            if (file_len >= 8)
            {
                uint8_t *dptr = reinterpret_cast<uint8_t*>(data.data());
                input.read(reinterpret_cast<char*>(dptr), file_len);
            }
            return data;
        }

        /// \brief read_ifd_offset  Reads offset to first IFD
        /// \param  input           Vector containing at least 8 bytes
        /// \returns                Offset to main IFD, or zero if format is not recognized
        uint32_t read_ifd_offset(std::vector<uint8_t> &data)
        {
            auto sz = data.size();
            if (sz < 8)
                return 0;
            auto dptr = reinterpret_cast<uint8_t *>(data.data());
            auto header = get_short(dptr, false);
            if (header == 0x4d4d && get_short(dptr + 2, true) == 42)
            {
                // Motorola format (network byte order or big endian) detected
                is_network_order = true;
            }
            else if (header == 0x4949 && get_short(dptr + 2, false) == 42)
            {
                // Intel format (little endian) detected
                is_network_order = false;
            }
            else
            {
                return 0;
            }

            return get_int(dptr + 4, is_network_order);
        }

        /// \brief read_ifd    Reads a tiff directory -- primary, secondary etc, or exif data
        /// \param data        Input vector
        /// \param offset      Starting offset to read the data from
        /// \returns           tuple of tags and offset to next IFD entry
        std::map<tiff_tags, tiff_tag_base>
            read_ifd(std::vector<uint8_t> &data, uint32_t &offset)
        {
            uint32_t data_len = (uint32_t)data.size();
            uint8_t *dptr = data.data();
            std::map<tiff_tags, tiff_tag_base> result;

            bool offset_is_valid = false;
            if (offset > 0 && offset + 2 <= data_len)
            {
                auto ifd_count = get_short(dptr + offset, is_network_order);
                offset += 2;
                // Check that data stream contains enough entries and offset to next directory
                if (offset + ifd_count * tiff_tag_size + 4 <= data_len)
                {
                    /// Read all entries
                    for (auto entry = 0; entry < ifd_count; entry++)
                    {
                        auto tag = tiff_tag_base(data, offset, is_network_order);
                        result[tag.tag()] = tag;
                        offset += tiff_tag_size;
                    }
                    // Read the link to next IFD
                    offset = get_int(dptr + offset, is_network_order);
                    offset_is_valid = true;
                }
            }
            offset = offset_is_valid ? offset : 0;
            return result;
        }

        /// Internal -- extract either the strip or tile data to internal arrays
        bool read_image_data(std::vector<uint8_t> &data, tiff_tags count_tag, tiff_tags offset_tag)
        {
            bool result = false;
            auto counts = tag_as_vector(count_tag);
            auto offsets = tag_as_vector(offset_tag);
            if (counts.size() != offsets.size())
                return result;

            auto dptr = data.data();
            uint32_t data_len = (uint32_t)data.size();

            for (size_t i = 0; i < counts.size(); i++)
            {
                uint32_t offset = offsets[i];
                uint32_t size = counts[i];
                if (offset + size < data_len)
                    _strips.push_back(std::vector<uint8_t>(dptr + offset, dptr + offset + size));
                result = true;
            }
            return result;
        }

        // Locates place in the tiff file for the image strip / tile data
        // Converts the strips/tiles as proper tags (currently only strip data is supported)
        void tiff_make_image_data(uint32_t &offset)
        {
            std::vector<uint32_t> offsets;
            std::vector<uint32_t> sizes;
            for (size_t i = 0; i < _strips.size(); i++)
            {
                auto &data = _strips[i];
                uint32_t strip_len = (uint32_t)data.size();
                offsets.push_back(offset);
                sizes.push_back(strip_len);
                offset += ((strip_len + 1) >> 1) << 1;
            }
            // Add the data either as strips or tiles
            // Behavior is unspecified, if both kind of image tags are present
            bool is_tiled = tag_as_vector(tag_tile_byte_counts).size() > 0;
            add_tag({
                tiff_tag<LONG>(is_tiled ? tag_tile_offsets : tag_strip_offsets, offsets),
                tiff_tag<LONG>(is_tiled ? tag_tile_byte_counts : tag_strip_byte_counts, sizes)
            });
        }

        /// Adds large data after current offset in the tiff file
        static void fix_tag_offsets(std::map<tiff_tags, tiff_tag_base> &tags, uint32_t &offset)
        {
            for (auto &tag : tags)
                tag.second.fix_offset(offset);
        }

        /// Adds exif data after current offset
        void rebuild_exif_data(uint32_t &offset)
        {
            uint16_t exif_size = (uint16_t)_exif.size();
            if (exif_size == 0)
            {
                erase_tag(tag_exif_ifd);
                return;
            }

            fix_tag_offsets(_exif, offset);

            add_tag(tiff_tag<LONG>(tag_exif_ifd, offset));
            offset += 6 + tiff_tag_size * exif_size;
        }

        /// Internal function testing if horizontal or vertical resampling factors are valid
        static uint16_t check_factor(uint32_t luma, uint32_t chroma)
        {
            if (luma == chroma)
                return 1;
            // Check it both ways to ensure no overflow happens
            if ((luma / 2 == chroma) && (luma % 2) == 0)
                return 2;
            if ((luma / 4 == chroma) && (luma % 4) == 0)
                return 4;
            throw std::invalid_argument("Chroma subsampling must be 1, 2 or 4");
        }

        /// Internal function to serialize and shift any image type
        /// to a contiguous buffer (of type uint8_t, uint16_t and uint32_t)
        /// with image shifted to maximum dynamic range
        template <typename T, typename U>
        void serialize_image(image<U> &input, int target_bpp)
        {
            // Allocate one strip worth of data
            _strips.clear();
            _strips.push_back(std::vector<uint8_t>(input._height * input._width * sizeof(T)));
            // Then reinterpret the buffer as typed image doing the conversion
            image<T> data(input._height, input._width, (T*)_strips[0].data());
            const int shift = sizeof(T) * 8 - target_bpp;
            data.foreach([shift](T &dst, U &src)
            {
                dst = (T)((T)src << shift);
            }, input);
        }

        /// Internal function to interleave, serialize and shift any image type
        /// to a contiguous buffer (of type uint8_t, uint16_t and uint32_t)
        /// with image shifted to maximum dynamic range
        template <typename T, typename U>
        void serialize_image(image<U>(&input)[3], int target_bpp, int hs = 1, int vs = 1)
        {
            // Scans vs*hs pixels from plane 0 for every pixel in planes 1 and 2
            int height = input[0]._height / vs;
            int width = input[0]._width / hs;
            int channels = hs*vs + 2;
            const int shift = sizeof(T) * 8 - target_bpp;

            // Allocate one strip worth of data -- then reinterpret it as typed image
            _strips.clear();
            _strips.push_back(std::vector<uint8_t>(height * width * channels * sizeof(T)));

            image<T> data(height, width * channels, (T*)_strips[0].data());

            for (int ch = 0; ch < channels; ch++)
            {
                auto input_channel = (ch < hs*vs) ?
                    input[0].subview(vs, hs, ch / vs, ch % vs) : input[ch - hs*vs + 1];
                auto output_channel = data.subview(1, channels, 0, ch);

                output_channel.foreach([shift](T &dst, U &src)
                {
                    dst = (T)((T)src << shift);
                }, input_channel);
            }
        }

        /// \brief  write_buffer    Writes image with pre-formatted data
        /// \param  out_stream      Output device/stream for the tiff file
        /// \param  width           Width of image
        /// \param  height          Height of image
        /// \param  bpp             Bits per pixel (8,16 or 32)
        /// \param  image_type      RGB, YUV or GRAY for bayer
        void write_buffer(
            std::ostream &out_stream,
            const int width, const int height,
            const int bpp, const int image_type)
        {
            switch (image_type)
            {
            case tag_photometric_interpretation_bayer:
                // 1 plane, bits per sample and sample format are scalars
                add_tag({
                    tiff_tag<SHORT>(tag_sample_format, tag_sample_format_unsigned),
                    tiff_tag<SHORT>(tag_bits_per_sample, (uint16_t)bpp),
                    tiff_tag<SHORT>(tag_samples_per_pixel, 1)
                });
                break;
            case tag_photometric_interpretation_yuv:
                // YUV -specific requirements -- tags 529-532
                add_tag({
                    tiff_tag<RATIONAL>(tag_reference_blackwhite,
                        { 0, 1, 255, 1, 128, 1, 255, 1, 128, 1, 255, 1 }),
                    tiff_tag<RATIONAL>(tag_ycbcr_coefficients, {299, 1000, 587, 1000, 114, 1000}),
                    tiff_tag<SHORT>(tag_ycbcr_positioning, tag_ycbcr_positioning_center)
                });
                // And if subsampling has not been given, add the default value
                if (tag_as_string(tag_ycbcr_subsampling) == "")
                    add_tag(tiff_tag<SHORT>(tag_ycbcr_subsampling, { 1, 1 }));

                // Fall through -- YUV and RGB share the planes == 3 related tags
            case tag_photometric_interpretation_rgb:
                // 3 planes, bits per sample and sample format are vectors
                add_tag({
                    tiff_tag<SHORT>(tag_sample_format, { 1, 1, 1 }),
                    tiff_tag<SHORT>(tag_bits_per_sample, std::vector<uint16_t>(3, (uint16_t)bpp)),
                    tiff_tag<SHORT>(tag_samples_per_pixel, 3)
                });
                break;
            default:
                throw std::invalid_argument("image_type must be RGB, YUV or Bayer");
            }

            add_tag({
                tiff_tag<SHORT>(tag_image_width, (uint16_t)width),
                tiff_tag<SHORT>(tag_image_height, (uint16_t)height),
                tiff_tag<SHORT>(tag_compression_method, tag_compression_method_uncompressed),
                tiff_tag<SHORT>(tag_photometric_interpretation, (uint16_t)image_type),
                tiff_tag<SHORT>(tag_orientation, tag_orientation_top_left),
                tiff_tag<SHORT>(tag_rows_per_strip, (uint16_t)height),
                tiff_tag<SHORT>(tag_planar_configuration, tag_planar_configuration_chunky)
            });

            write(out_stream);
        }

        // Write 8-byte Tiff header with offset to first IFD directory
        void write_tiff_header(std::ostream &out_stream, uint32_t offset)
        {
            put_short(out_stream, is_network_order ? 0x4D4D : 0x4949, true);
            put_short(out_stream, 42, is_network_order);
            put_int(out_stream, offset, is_network_order);
        };

        /// Write the serialized image buffer and an additional padding
        /// to next uint16_t word boundary -- TODO: fix byte order (for 16-bit bigendian)
        void write_image_strips(std::ostream &out_stream)
        {
            for (size_t i = 0; i < _strips.size(); i++)
            {
                auto &data = _strips[i];
                uint32_t strip_len = (uint32_t)data.size();
                out_stream.write((char *)data.data(), strip_len);
                if (strip_len & 1)
                    out_stream.put(0);
            }
        }

        /// Serializes all additional data from tags to file, that didn't fit in the 4-byte offset
        /// All offsets should be even, so the data will be written in chunks of 16-bits
        void write_additional_data(std::ostream &stream, std::map<tiff_tags, tiff_tag_base> &dir)
        {
            for (auto &x : dir)
            {
                auto tag = x.second;
                tag.write_additional_data(stream, is_network_order);
            }
        }

        /// Write number of tags, full ifd tag dictionary and offset the next directory (0)
        ///         +-----+--//---+-------+-----+----------------+
        /// IFD #0: |count| count x IFD_ENTRIES |Offset to IFD #1|
        ///         +-----+--//---+-------+-----+----------------+
        void write_ifd(std::ostream &out_stream, std::map<tiff_tags, tiff_tag_base> &dir)
        {
            put_short(out_stream, (uint16_t)dir.size(), is_network_order);
            for (auto &x : dir)
            {
                auto tag = x.second;
                tag.write_pod(out_stream, is_network_order);
            }
            put_int(out_stream, 0, is_network_order);
        }
    };
};