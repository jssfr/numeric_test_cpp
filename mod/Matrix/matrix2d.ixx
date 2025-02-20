module;

#include <atomic>
#include <stdexcept>
#include<format>
#include<fmt/core.h>
#include "std_input.hpp"
#include "output.hpp"
#include<ranges>

export module Matrix:matrix2d;

import Config;
import math;

namespace jf::matrix {

    // @link https://www.youtube.com/@HomoSiliconiens
template <typename ElementType>
using aligned_safe_vector_t =
    oneapi::tbb::concurrent_vector<ElementType, oneapi::tbb::cache_aligned_allocator<ElementType>>;

template <typename ElementType>
using aligned_fast_vector_t =
    std::vector<ElementType, oneapi::tbb::cache_aligned_allocator<ElementType>>;

template <typename ElementType>
using scalable_safe_vector_t =
    oneapi::tbb::concurrent_vector<ElementType, oneapi::tbb::scalable_allocator<ElementType>>;

template <typename ElementType>
using scalable_fast_vector_t =
    std::vector<ElementType, oneapi::tbb::scalable_allocator<ElementType>>;

template <typename ElementType, template <typename...> class ContainerType = std::vector,
          template <typename> class AllocatorType = std::allocator>
class matrix_2d {
   public:
    using value_type = ElementType;
    using allocator_t = AllocatorType<value_type>;
    using container_t = ContainerType<value_type, allocator_t>;
    // unidimentional vector is faster then vector of vectors  aka vector<vector<> >
   private:
    size_t m_rows, m_cols;
    container_t m_array;

    void invalidate() {
        this->m_rows = 0;
        this->m_cols = 0;
    }

    template <typename IndexType>
    inline bool row_range_valid(IndexType row) const {
        return (row >= 0 && row < this->m_rows);
    }

    template <typename IndexType>
    inline bool column_range_valid(IndexType col) const {
        return (col >= 0 && col < this->m_cols);
    }
    template <typename IndexType>
    inline bool is_range_valid(IndexType row, IndexType col) const noexcept {
        return row_range_valid(row) && column_range_valid(col);
    }

    template <typename IndexElementType, template <typename, typename...> class IndexContainerType,
              typename... Types>
    ElementType term(const IndexContainerType<IndexElementType, Types...>& index) const {
        ElementType rlt = 1.0;

        for (size_t i = 0; i < this->m_rows; ++i) rlt *= this->operator()(i, (size_t)index[i]);

        //return jf::ncrnpr::sgn(index) == 1? rlt : -rlt;
        return rlt;
    }

   public:
    ElementType determinant_leibniz_serial() const {
        if (this->m_rows == 0)
            return ElementType{};
        else if (this->m_rows == 1)
            return this->m_array[0];
        else if (this->m_rows == 2) {
            return this->operator()(0, 0) * this->operator()(1, 1) -
                   this->operator()(0, 1) * this->operator()(1, 0);
        } else {
            ElementType det = ElementType{};

            std::vector<size_t> index;
            std::generate_n(std::back_inserter(index), this->m_rows,
                            [count = 0]() mutable { return count++; });

            do {
                if (jf::ncrnpr::sgn(index) > 0)
                    det += this->term(index);
                else
                    det -= this->term(index);
            } while (std::next_permutation(index.begin(), index.end()));

            return det;
        }
    }

    ElementType determinant_leibniz_parallel() const {

        size_t n = this->m_rows, r;

        if (n < 6)
            r = 1;
        else if (n < 10)
            r = 2;
        else
            r = 3;

        auto max_permu = jf::ncrnpr::nPr(n, r);

        std::atomic<ElementType> det{};  // share among multiple threads;

        auto inner_loop = [&](auto mth) {
            auto index = jf::ncrnpr::enum_permu_remainder(n, r, mth);

            ElementType partial_det = ElementType{};

            do {
                if (jf::ncrnpr::sgn(index) == 1)
                    partial_det += this->term(index);


                else
                    partial_det -= this->term(index);
            } while (std::next_permutation(index.begin() + r, index.end()));

            // det+= partial_det; // works but doesnt work ! we working with mult threads so looks like working but doenst

            ElementType old_det = det, new_det = old_det + partial_det;

            while (!det.compare_exchange_strong(old_det, new_det)) {
                old_det = det, new_det = old_det + partial_det;
            }
        };

        oneapi::tbb::parallel_for(size_t{}, (size_t)max_permu, inner_loop);

        return det;
    }
    template <typename IndexType>
    matrix_2d minor(IndexType row, IndexType column) const {
        matrix_2d mm{this->m_rows - 1, this->m_cols - 1};
        size_t ii, jj;

        for (size_t i = 0; i < this->m_rows; ++i) {
            if (i == row) continue;

            ii = (i < row) ? i : i - 1;

            for (size_t j = 0; j < this->m_cols; ++j) {
                if (j == column) continue;
                jj = (j < column) ? j : j - 1;

                mm(ii, jj) = this->operator()(i, j);
            }
        }

        return mm;
    }

    ElementType determinant_laplace_serial() const {
        if (this->m_rows == 0)
            return {};
        else if (this->m_rows == 1)
            return this->m_array[0];
        else if (this->m_rows == 2) {
            return this->operator()(0, 0) * this->operator()(1, 1) -
                   this->operator()(0, 1) * this->operator()(1, 0);
        } else {
            ElementType det = ElementType{};

            for (size_t j = 0; j < this->m_cols; ++j) {
                auto mm = this->minor(size_t{}, j);
                if (j % 2)  // odd
                    det -= this->operator()(size_t{}, j) *
                           mm.determinant_laplace_serial();  // co-factor
                else                                         // even
                    det += this->operator()(size_t{}, j) *
                           mm.determinant_laplace_serial();  // co-factor
            }  // mm.determinant_laplace_serial() == Mathematical Induction,
            // Recurrency Relation,
            // Recursion or basicly Reduction; (n) -> (n-1, n-2, ...)
            return det;
        }
    }

    ElementType determinant_laplace_parallel_atomic() const {
        if (this->m_rows == 0)
            return {};
        else if (this->m_rows == 1)
            return this->m_array[0];
        else if (this->m_rows == 2) {
            return this->operator()(0, 0) * this->operator()(1, 1) -
                   this->operator()(0, 1) * this->operator()(1, 0);
        } else {
            std::atomic<ElementType> det{};

            auto handle = [&](auto& range) {  // auto& ,aps to oneapi::tbb::blocked_range<size_t>&
                for (auto j = range.begin(); j < range.end(); ++j) {
                    auto mm = this->minor(size_t{}, j);
                    auto cofactor = this->operator()(size_t{}, j) *
                                    mm.determinant_laplace_serial();  // co-factor

                    ElementType old_det = det;
                    ElementType new_det = (j % 2) ? old_det - cofactor : old_det + cofactor;
                    while (!det.compare_exchange_strong(old_det, new_det)) {
                        old_det = det;
                        new_det = (j % 2) ? old_det - cofactor : old_det + cofactor;
                    }
                }

            };

            oneapi::tbb::parallel_for(oneapi::tbb::blocked_range{size_t{}, this->m_cols}, handle);

            return det;
        }
    }

    ElementType determinant_laplace_parallel_reduce() const {
        if (this->m_rows == 0)
            return {};
        else if (this->m_rows == 1)
            return this->m_array[0];
        else if (this->m_rows == 2) {
            return this->operator()(0, 0) * this->operator()(1, 1) -
                   this->operator()(0, 1) * this->operator()(1, 0);
        } else {
            auto handle = [&](auto& range, auto det) {
                for (auto j = range.begin(); j != range.end(); ++j) {
                    auto mm = this->minor(size_t{}, j);
                    auto cofactor = this->operator()(size_t{}, j) *
                                    mm.determinant_laplace_serial();  // co-factor

                    det = (j % 2) ? (det - cofactor) : (det + cofactor);
                }
                return det;
            };

            auto sum_up = [](auto left_det, auto right_det) { return left_det + right_det; };

            return oneapi::tbb::parallel_reduce(oneapi::tbb::blocked_range{size_t{}, this->m_cols},
                                                ElementType{}, handle, sum_up);
        }
    }

    size_t index_base() const noexcept { return 0; }

    inline bool empty() const noexcept { return this->m_rows == 0; }

    inline bool not_empty() const noexcept { return this->m_rows != 0; }

    inline size_t rows() const noexcept { return this->m_rows; }

    inline size_t columns() const noexcept { return this->m_cols; }

    matrix_2d() noexcept : m_rows{}, m_cols{}, m_array{} {}

    template <typename SizeType>
    matrix_2d(SizeType rows, SizeType cols)
        : m_rows{static_cast<size_t>(rows)},
          m_cols{static_cast<size_t>(cols)},
          m_array(static_cast<size_t>(rows * cols)) {
              std::ranges::for_each(this->m_array, [](auto& i){ i = 0; });
          }

    matrix_2d(const matrix_2d&) = default;
    matrix_2d& operator=(const matrix_2d&) = default;

    matrix_2d(matrix_2d&& rhs) noexcept
        : m_rows{rhs.m_rows}, m_cols{rhs.m_cols}, m_array{std::move(rhs.m_array)} {
        rhs.invalidate();
    }

    matrix_2d& operator=(matrix_2d&& rhs) noexcept {
        if (this != std::addressof(rhs) && rhs.not_empty()) {
            this->m_rows = rhs.m_rows;
            this->m_cols = rhs.m_cols;
            this->m_array = std::move(rhs.m_array);
            rhs.invalidate();
        }
        return *this;
    }

    value_type& operator()(size_t row, size_t col) {
        if (row_range_valid(row) && column_range_valid(col)) {
            return this->m_array[row * m_cols + col];
        } else {
            std::cerr << "File name [" << __FILE__ << "]\n"
                      << "Line no: " << __LINE__ << ", index out of range\n"
                      << "Valid range: 0 <= row < " << this->m_rows << ", requested row: " << row
                      << "\n"
                      << "Valid range: 0 <= column < " << this->m_cols
                      << ", requested column: " << col << "\n";
            std::abort();
        }
    }

    const value_type& operator()(size_t row, size_t col) const {
        if (row_range_valid(row) && column_range_valid(col)) {
            return this->m_array[row * m_cols + col];
        } else {
            std::cerr << "File name [" << __FILE__ << "]\n"
                      << "Line no: " << __LINE__ << ", index out of range\n"
                      << "Valid range: 0 <= row < " << this->m_rows << ", requested row: " << row
                      << "\n"
                      << "Valid range: 0 <= column < " << this->m_cols
                      << ", requested column: " << col << "\n";
            std::abort();
        }
    }

    void set_value(std::vector<value_type> cntr){
        if(cntr.empty()) return;

        // size_t x{1};
        // size_t y{1};
        // if(cntr.size() >= this->m_rows * this->m_cols){
        //     x = this->m_rows;
        //     y = this->m_cols;
        // }else if(cntr <= this->m_cols){
        //     y = cntr.size();
        // }else{
        //     y = this->m_cols;  
        //     x = cntr.size() - y; // 
        // }
        // for(std::size_t i{}; i < x; ++i){
        //     for(std::size_t j{}; j < y; ++j){
        //         this->operator()(i, j) = v1++;
        //     }
        // }
        if(cntr.size() >= this->m_rows * this->m_cols){
            for(size_t i{}; i < this->m_array.size(); ++i){
                this->m_array[i] = cntr[i];
            }
        }else{
            for(size_t i{}; i < cntr.size(); ++i){
                this->m_array[i] = cntr[i];
            }    
        }

    }

    matrix_2d& operator*=(value_type scalar) {
        oneapi::tbb::parallel_for_each(m_array.begin(), m_array.end(),
                                       [scalar](auto& ele) { ele *= scalar; });

        return *this;
    }


    matrix_2d& operator/=(value_type scalar) {
        if (scalar == 0) { throw std::invalid_argument("Cannot divide by 0."); }

        *this *= value_type{1} / scalar;

        return *this;
    }

    friend matrix_2d operator*(const matrix_2d& m, value_type scalar) {
        auto mm = m;
        mm *= scalar;
        return mm;
    }

    friend matrix_2d operator*(value_type scalar, const matrix_2d& m) {
        auto mm = m;
        mm *= scalar;
        return mm;
    }

    friend matrix_2d&& operator*(matrix_2d&& m, value_type scalar) {
        m *= scalar;
        return std::move(m);
    }

    friend matrix_2d&& operator*(value_type scalar, matrix_2d&& m) {
        m *= scalar;
        return std::move(m);
    }

    friend matrix_2d operator*( const matrix_2d& m1, const matrix_2d& m2) {
        
        if(!(m1.columns() == m2.rows())){ throw std::invalid_argument("Cannot multiply matrices."); }
/*
        | a v c|        |q|      |p|
        | d e f|    *   |w|    = |l|
        | g h i|3x3     |j|3x1   |m|3x1
*/
        matrix_2d mm{m1.rows(), m2.columns()};

        oneapi::tbb::parallel_for(size_t{}, m1.rows(), [&](auto i){
             for(size_t j{}; j < mm.columns(); ++j){
                for (size_t k{}; k < mm.rows(); ++k) {
                 mm(i, j) += m1(i, k) * m2(k, j);
                 }
            }
        });

        return mm;
    }
    friend matrix_2d operator/(const matrix_2d& m, value_type scalar) {
        if (scalar == 0) { throw std::invalid_argument("Cannot divide by 0."); }

        auto mm = m;
        mm /= scalar;
        return mm;
    }

    friend matrix_2d&& operator/(matrix_2d&& m, value_type scalar) {
        if (scalar == 0) { throw std::invalid_argument("Cannot divide by 0."); }

        m /= scalar;
        return std::move(m);
    }
};

export template <typename ElementType>
using aligned_safe_matrix_2d_t =
    matrix_2d<ElementType, oneapi::tbb::concurrent_vector, oneapi::tbb::cache_aligned_allocator>;

export template <typename ElementType>
using aligned_fast_matrix_2d_t =
    matrix_2d<ElementType, std::vector, oneapi::tbb::cache_aligned_allocator>;

export template <typename ElementType>
using scalable_safe_matrix_2d_t =
    matrix_2d<ElementType, oneapi::tbb::concurrent_vector, oneapi::tbb::scalable_allocator>;

export template <typename ElementType>
using scalable_fast_matrix_2d_t =
    matrix_2d<ElementType, std::vector, oneapi::tbb::scalable_allocator>;

export template <typename ElementType>
using mat =
    matrix_2d<ElementType, std::vector, oneapi::tbb::scalable_allocator>;

export using fmat = matrix_2d<float, std::vector, oneapi::tbb::scalable_allocator>;

export using dmat = matrix_2d<double, std::vector, oneapi::tbb::scalable_allocator>;

export using imat = matrix_2d<int, std::vector, oneapi::tbb::scalable_allocator>;

} //namespace jf::matrix
 
 
// Função de formatação para fmt
template <typename ElementType, template <typename...> class ContainerType,
          template <typename> class AllocatorType>
struct fmt::formatter<jf::matrix::matrix_2d<ElementType, ContainerType, AllocatorType>> {
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const jf::matrix::matrix_2d<ElementType, ContainerType, AllocatorType>& m,
                FormatContext& ctx) const{
        auto out = ctx.out();
        if (m.empty()) { return fmt::format_to(out, "[ ]"); }

        fmt::format_to(out, "\n");
        for (size_t i = m.index_base(); i < m.rows(); ++i) {
            fmt::format_to(out, "[");
            for (size_t j = m.index_base(); j < m.columns(); ++j) {
                fmt::format_to(out, "{}", m(i, j));
                if (j + 1 < m.columns()) {
                    fmt::format_to(out, ", ");
                } else {
                    fmt::format_to(out, "]");
                }
            }
            fmt::format_to(out, "\n");
        }
        return out;
    }
};

template <typename ElementType, template <typename...> class ContainerType,
          template <typename> class AllocatorType>
struct std::formatter<jf::matrix::matrix_2d<ElementType, ContainerType, AllocatorType>> {
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const jf::matrix::matrix_2d<ElementType, ContainerType, AllocatorType>& m,
                FormatContext& ctx) const{
        auto out = ctx.out();
        if (m.empty()) { return std::format_to(out, "[ ]"); }

        std::format_to(out, "\n");
        for (size_t i = m.index_base(); i < m.rows(); ++i) {
            std::format_to(out, "[");
            for (size_t j = m.index_base(); j < m.columns(); ++j) {
                std::format_to(out, "{}", m(i, j));
                if (j + 1 < m.columns()) {
                    std::format_to(out, ", ");
                } else {
                    std::format_to(out, "]");
                }
            }
            std::format_to(out, "\n");
        }
        return out;
    }
};

