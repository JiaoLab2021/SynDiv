#ifndef CHUNKREADER_HPP
#define CHUNKREADER_HPP
#include <fstream>
#include <string.h>
#include "get_time.hpp"

class ChunkReader {
public:
    ChunkReader(const std::string& file_path, const size_t chunk_size = 1024 * 1024)  // 默认1Mb的缓存
        : file_path_(file_path), chunk_size_(chunk_size), chunk_buffer_(new char[chunk_size_]), buffer_pos_(0), buffer_end_(0) {
            file_.open(file_path_, std::ifstream::binary);
            if (!file_.is_open()) {
                std::cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "'" << file_path_ << "': No such file or directory or possibly reached the maximum open file limit. You can set 'ulimit -n' to a larger value to continue." << endl;
                exit(1);
            }
    }

    ~ChunkReader() {
        if (file_.is_open()) {
            file_.close();
        }
    }

    bool read_line(std::string& line) {
        line.clear();

        while (true) {
            // Fill the buffer if it's empty
            if (buffer_pos_ == buffer_end_) {
                if (!read_chunk()) {
                    // 如果到达文件结尾，但是没有换行符，返回最后一部分内容。不然少读一行
                    if (!line.empty())
                    {
                        return true;
                    }
                    return false;
                }
            }
            
            // Extract a line from the buffer
            const char* line_begin = buffer_ + buffer_pos_;
            const char* line_end = static_cast<const char*>(memchr(line_begin, '\n', buffer_end_ - buffer_pos_));
            if (line_end) {
                line.append(line_begin, line_end - line_begin);
                buffer_pos_ += static_cast<size_t>(line_end - line_begin) + 1;
                return true;
            } else {
                line.append(line_begin, buffer_end_ - buffer_pos_);
                buffer_pos_ = buffer_end_;
            }
        }

        return false;
    }

private:
    bool read_chunk() {
        if (!file_.good()) {
            return false;
        }
        file_.read(chunk_buffer_.get(), chunk_size_);
        buffer_pos_ = 0;
        buffer_end_ = file_.gcount();
        if (buffer_end_ == 0) {
            return false;
        }
        buffer_ = chunk_buffer_.get();
        return true;
    }

private:
    std::string file_path_;
    size_t chunk_size_;
    std::unique_ptr<char[]> chunk_buffer_;
    char* buffer_;
    size_t buffer_pos_;
    size_t buffer_end_;
    std::ifstream file_;
};

#endif