module;
#ifndef USING_IMPORT_STD_MOD
#  include <array>
#  include <atomic>
#  include <format>
#  include <iostream>
#  include <ranges>
#  include <stdexcept>
#  include <type_traits>
// #  include "output.hpp"
#  include "std_input.hpp"
#  ifdef __clang__
#    ifdef USING_TBBLIB
#      include <oneapi/tbb.h>
#    endif
#  endif
#else
//
#  ifdef __clang__
#    ifdef USING_TBBLIB
#      include <oneapi/tbb.h>
#    endif
#  endif
#endif

export module math:matrix2d;

#ifdef USING_IMPORT_STD_MOD
import std;
#endif
import Config;
import parallel;
import :ncrnpr;

namespace jf::matrix
{
#ifdef __clang__
#  ifdef USING_TBBLIB
namespace parallel = oneapi::tbb;  // tbb is generating TU-local problems with gcc
#  else
namespace parallel = jf::par;
#  endif
#else
namespace parallel = jf::par;
#endif

/// @link https://www.youtube.com/@HomoSiliconiens
template<typename ElementType, std::size_t Rows, std::size_t Cols,
    template<typename...> class ContainerType = std::vector,
    template<typename> class AllocatorType = std::allocator>
  requires(Rows > 0 && Cols > 0)
class matrix_2d
{
public:
  using value_type = ElementType;
  using allocator_t = AllocatorType<value_type>;
  using container_t = ContainerType<value_type, allocator_t>;

private:
  std::size_t m_rows, m_cols;
  container_t m_array;

  void invalidate()
  {
    this->m_rows = 0;
    this->m_cols = 0;
  }

  template<typename IndexType> inline bool row_range_valid(IndexType row) const
  {
    return (row >= 0 && row < this->m_rows);
  }

  template<typename IndexType> inline bool column_range_valid(IndexType col) const
  {
    return (col >= 0 && col < this->m_cols);
  }

  template<typename IndexType>
  inline bool is_range_valid(IndexType row, IndexType col) const noexcept
  {
    return row_range_valid(row) && column_range_valid(col);
  }

  template<typename IndexElementType, template<typename, typename...> class IndexContainerType,
      typename... Types>
  ElementType term(const IndexContainerType<IndexElementType, Types...>& index) const
  {
    auto rlt = static_cast<ElementType>(1);

    for (std::size_t i {}; i < this->m_rows; ++i)
      rlt *= this->operator()(i, (std::size_t)index[i]);

    return rlt;
  }

public:
  ElementType determinant_leibniz_serial() const
  {
    if constexpr (this->m_rows == 0)
      return ElementType {};
    else if constexpr (this->m_rows == 1)
      return this->m_array[0];
    else if constexpr (this->m_rows == 2) {
        return this->operator()(0, 0) * this->operator()(1, 1)
            - this->operator()(0, 1) * this->operator()(1, 0);
    } else {
        ElementType det = ElementType {};

        // std::vector<std::size_t> index;
        // std::generate_n(std::back_inserter(index), this->m_rows,
        //     [count = 0]() mutable { return count++; });
        std::vector<std::size_t> index { std::views::iota(0ul) | std::views::take(this->m_rows)
          | std::ranges::to<std::vector>() };

        do {
            if (jf::ncrnpr::sgn(index) > 0)
              det += this->term(index);
            else
              det -= this->term(index);
        } while (std::next_permutation(index.begin(), index.end()));

        return det;
      }
  }

  ElementType determinant_leibniz_parallel() const
  {
    std::size_t n = this->m_rows, r;

    if (n < 6)
      r = 1;
    else if (n < 10)
      r = 2;
    else
      r = 3;

    auto max_permu = jf::ncrnpr::nPr(n, r);

    std::atomic<ElementType> det {};  // share among multiple threads;

    auto inner_loop = [&](auto mth)
    {
      auto index = jf::ncrnpr::enum_permu_remainder(n, r, mth);

      ElementType partial_det = ElementType {};

      do {
          if (jf::ncrnpr::sgn(index) == 1)
            partial_det += this->term(index);
          else
            partial_det -= this->term(index);
      } while (std::next_permutation(index.begin() + r, index.end()));

      ElementType old_det = det, new_det = old_det + partial_det;

      while (!det.compare_exchange_strong(old_det, new_det)) {
          old_det = det, new_det = old_det + partial_det;
        }
    };

    parallel::parallel_for(std::size_t {}, static_cast<std::size_t>(max_permu), inner_loop);

    return det;
  }

  template<typename IndexType> auto minor(IndexType row, IndexType column) const
  {
    const std::size_t rows = Rows - 1;
    const std::size_t cols = Cols - 1;
    matrix_2d<std::size_t, rows, cols, ContainerType, AllocatorType> mm;
    std::size_t ii, jj;

    for (std::size_t i = 0; i < this->m_rows; ++i) {
        if (i == row) continue;

        ii = (i < row) ? i : i - 1;

        for (std::size_t j = 0; j < this->m_cols; ++j) {
            if (j == column) continue;
            jj = (j < column) ? j : j - 1;

            mm(ii, jj) = this->operator()(i, j);
          }
      }

    return mm;
  }

  ElementType determinant_laplace_serial() const
  {
    if constexpr (Rows == 0)
      return {};
    else if constexpr (Rows == 1)
      return this->m_array[0];
    else if constexpr (Rows == 2) {
        return this->operator()(0, 0) * this->operator()(1, 1)
            - this->operator()(0, 1) * this->operator()(1, 0);
    } else {
        ElementType det = ElementType {};

        for (std::size_t j = 0; j < this->m_cols; ++j) {
            auto mm = this->minor(std::size_t {}, j);
            if (j % 2)  // odd
              det -= this->operator()(std::size_t {}, j)
                  * mm.determinant_laplace_serial();  // co-factor
            else  // even
              det += this->operator()(std::size_t {}, j)
                  * mm.determinant_laplace_serial();  // co-factor
          }  // mm.determinant_laplace_serial() == Mathematical Induction,
        // Recurrency Relation,
        // Recursion or basicly Reduction; (n) -> (n-1, n-2, ...)
        return det;
      }
  }

  ElementType determinant_laplace_parallel_atomic() const
  {
    if (this->m_rows == 0) return {};
    if (this->m_rows == 1) return this->m_array[0];
    if (this->m_rows == 2) {
        return this->operator()(0, 0) * this->operator()(1, 1)
            - this->operator()(0, 1) * this->operator()(1, 0);
    }

    std::atomic<ElementType> det {};

    auto handle = [&](auto &range) 
    {                     // auto& ,aps to
                         // oneapi::tbb::blocked_range<std::size_t>& or
                         // jf::par::blocked_range<std::size_t>&
      for (auto j = range.begin(); j < range.end(); ++j) {
        auto mm = this->minor(std::size_t{}, j);
        auto cofactor = this->operator()(std::size_t{}, j) *
                          mm.determinant_laplace_serial(); // co-factor

        ElementType old_det = det;
        ElementType new_det = (j % 2) ? old_det - cofactor : old_det + cofactor;

        while (!det.compare_exchange_strong(old_det, new_det)) {
            old_det = det;
            new_det = (j % 2) ? old_det - cofactor : old_det + cofactor;
        }
      }

    };

    parallel::parallel_for(parallel::blocked_range { std::size_t {}, this->m_cols }, handle);

    return det;
  }

  ElementType determinant_laplace_parallel_reduce() const
  {
    if (this->m_rows == 0) return {};
    if (this->m_rows == 1) return this->m_array[0];
    if (this->m_rows == 2) {
        return this->operator()(0, 0) * this->operator()(1, 1)
            - this->operator()(0, 1) * this->operator()(1, 0);
    }

    auto handle = [&](auto& range, auto det)
    {
      for (auto j = range.begin(); j != range.end(); ++j) {
          auto mm = this->minor(std::size_t {}, j);
          auto cofactor = this->operator()(std::size_t {}, j)
              * mm.determinant_laplace_serial();  // co-factor

          det = (j % 2) ? (det - cofactor) : (det + cofactor);
        }
      return det;
    };

    auto sum_up = [](auto left_det, auto right_det) { return left_det + right_det; };

    return parallel::parallel_reduce(parallel::blocked_range { std::size_t {}, this->m_cols },
        ElementType {}, handle, sum_up);
  }

  inline bool empty() const noexcept { return this->m_rows == 0; }

  inline std::size_t rows() const noexcept { return this->m_rows; }

  inline std::size_t columns() const noexcept { return this->m_cols; }

  matrix_2d() noexcept
      : m_rows { Rows }
      , m_cols { Cols }
      , m_array(Rows * Cols)
  {
    std::ranges::fill(this->m_array, 0);
  }

  matrix_2d(std::array<std::array<ElementType, Cols>, Rows> cntr)
  {
    this->m_array = cntr | std::views::join | std::ranges::to<std::vector>();
  }

  matrix_2d(const matrix_2d&) = default;
  matrix_2d& operator=(const matrix_2d&) = default;

  matrix_2d(matrix_2d&& rhs) noexcept
      : m_rows { rhs.m_rows }
      , m_cols { rhs.m_cols }
      , m_array { std::move(rhs.m_array) }
  {
    rhs.invalidate();
  }

  matrix_2d& operator=(matrix_2d&& rhs) noexcept
  {
    if (this != std::addressof(rhs) && !rhs.empty()) {
        this->m_rows = rhs.m_rows;
        this->m_cols = rhs.m_cols;
        this->m_array = std::move(rhs.m_array);
        rhs.invalidate();
    }
    return *this;
  }

  value_type& operator()(std::size_t row, std::size_t col)
  {
    if (row_range_valid(row) && column_range_valid(col)) {
        return this->m_array[row * m_cols + col];
    } else {
        std::cerr << "File name [" << __FILE__ << "]\n"
                  << "Line no: " << __LINE__ << ", index out of range\n"
                  << "Valid range: 0 <= row < " << this->m_rows << ", requested row: " << row
                  << "\n"
                  << "Valid range: 0 <= column < " << this->m_cols << ", requested column: " << col
                  << "\n";
        std::abort();
      }
  }

  const value_type& operator()(std::size_t row, std::size_t col) const
  {
    if (row_range_valid(row) && column_range_valid(col)) {
        return this->m_array[row * m_cols + col];
    } else {
        std::cerr << "File name [" << __FILE__ << "]\n"
                  << "Line no: " << __LINE__ << ", index out of range\n"
                  << "Valid range: 0 <= row < " << this->m_rows << ", requested row: " << row
                  << "\n"
                  << "Valid range: 0 <= column < " << this->m_cols << ", requested column: " << col
                  << "\n";
        std::abort();
      }
  }

  void set_value(std::array<std::array<ElementType, Cols>, Rows> cntr)
  {
    this->m_array = cntr | std::views::join | std::ranges::to<std::vector>();
  }

  // void set_value(std::vector<value_type> cntr)
  // {
  //   if (cntr.empty()) return;
  //
  //   if (cntr.size() >= this->m_rows * this->m_cols) {
  //       for (std::size_t i {}; i < this->m_array.size(); ++i) {
  //           this->m_array[i] = cntr[i];
  //         }
  //   } else {
  //       for (std::size_t i {}; i < cntr.size(); ++i) {
  //           this->m_array[i] = cntr[i];
  //         }
  //     }
  // }

  template<typename S>
    requires std::is_arithmetic_v<S>
  matrix_2d& operator*=(S scalar)
  {
    parallel::parallel_for_each(m_array.begin(), m_array.end(),
        [scalar](auto& ele) { ele *= scalar; });

    return *this;
  }

  template<typename S>
    requires std::is_arithmetic_v<S>
  matrix_2d& operator/=(S scalar)
  {
    if (scalar == 0) {
        throw std::invalid_argument("Cannot divide by 0.");
    }

    *this *= value_type { 1 } / scalar;

    return *this;
  }

  /// houseHolder method for orthogonal and upper triangular matrices for QR algorithm
  /// @link https://en.wikipedia.org/wiki/Householder_transformation
  auto houseHolder() -> std::tuple<matrix_2d, matrix_2d>
  {
    matrix_2d R = *this;
    const std::size_t rows { Rows };
    const std::size_t cols { Cols };
    matrix_2d<ElementType, rows, cols, ContainerType, AllocatorType> Q;

    auto norm = [&](const matrix_2d& v, auto k)
    {
      value_type sum {};
      for (auto i = k; i < v.rows(); ++i) {
          sum += v(i, k) * v(i, k);
        }
      return std::sqrt(sum);
    };

    for (std::size_t i {}; i < rows; ++i) {
        Q(i, i) = static_cast<value_type>(1);
      }

    for (std::size_t k = 0; k < cols; ++k) {
        // Householder's vector
        std::vector<value_type> v(rows, static_cast<value_type>(0));

        value_type alpha = -std::copysign(norm(R, k), R(k, k));
        value_type r = std::sqrt(0.5 * (alpha * alpha - R(k, k) * alpha));

        if (r == static_cast<value_type>(0)) continue;

        v[k] = (R(k, k) - alpha) / (2 * r);
        for (std::size_t i = k + 1; i < rows; ++i) {
            v[i] = R(i, k) / (2 * r);
          }

        // update R
        for (std::size_t j = k; j < cols; ++j) {
            value_type dot {};
            for (std::size_t i = k; i < rows; ++i) {
                dot += v[i] * R(i, j);
              }
            for (std::size_t i = k; i < rows; ++i) {
                R(i, j) -= 2 * v[i] * dot;
              }
          }

        // update Q
        for (std::size_t i = 0; i < rows; ++i) {
            value_type dot {};
            for (std::size_t j = k; j < cols; ++j) {
                dot += v[j] * Q(i, j);
              }
            for (std::size_t j = k; j < cols; ++j) {
                Q(i, j) -= 2 * v[j] * dot;
              }
          }
      }

    return std::tuple { std::move(Q), std::move(R) };
  }

  auto houseHolder_parallel() -> std::tuple<matrix_2d, matrix_2d>
  {
    matrix_2d R = *this;
    matrix_2d Q { R.identity() };  // R.rows(), R.columns() };

    std::size_t rows { this->m_rows };
    std::size_t cols { this->m_cols };

    auto norm = [&](const matrix_2d& v, auto k)
    {
      value_type sum {};
      for (auto i = k; i < v.rows(); ++i) {
          sum += v(i, k) * v(i, k);
        }
      return std::sqrt(sum);
    };

    // for (std::size_t i {}; i < rows; ++i) {
    //     Q(i, i) = static_cast<value_type>(1);
    //   }

    auto sum_up = [](auto left_sum, auto right_sum) { return left_sum + right_sum; };

    for (std::size_t k = 0; k < cols; ++k) {
        // Householder's vector
        std::vector<value_type> v(rows, static_cast<value_type>(0));

        value_type alpha = -std::copysign(norm(R, k), R(k, k));
        value_type r = std::sqrt(0.5 * (alpha * alpha - R(k, k) * alpha));

        if (r == static_cast<value_type>(0)) continue;

        v[k] = (R(k, k) - alpha) / (2 * r);

        parallel::parallel_for(parallel::blocked_range { k + 1, rows },
            [&](auto& range)
            {
              for (std::size_t i = range.begin(); i < range.end(); ++i) {
                  v[i] = R(i, k) / (2 * r);
                }
            });

        // update R
        for (std::size_t j = k; j < cols; ++j) {
            value_type dot = parallel::parallel_reduce(
                parallel::blocked_range { k, rows }, value_type {},
                [&](auto& range, value_type sum)
                {
                  for (std::size_t i = range.begin(); i < range.end(); ++i)
                    sum += v[i] * R(i, j);
                  return sum;
                },
                sum_up);

            parallel::parallel_for(parallel::blocked_range { k, rows },
                [&](auto& range)
                {
                  for (std::size_t i = range.begin(); i < range.end(); ++i)
                    R(i, j) -= 2 * v[i] * dot;
                });
          }

        // update Q
        parallel::parallel_for(parallel::blocked_range { std::size_t {}, rows },
            [&](auto& range)
            {
              for (std::size_t i = range.begin(); i < range.end(); ++i) {
                  value_type dot {};
                  for (std::size_t j = k; j < cols; ++j) {
                      dot += v[j] * Q(i, j);
                    };

                  for (std::size_t j = k; j < cols; ++j) {
                      Q(i, j) -= 2 * v[j] * dot;
                    }
                }
            });
      }
    return std::tuple { std::move(Q), std::move(R) };
  }

  auto identity() -> matrix_2d
  {
    const std::size_t rows = Rows;
    const std::size_t cols = Cols;

    matrix_2d<value_type, rows, cols, ContainerType, AllocatorType> mat;
    for (std::size_t i {}; i < this->m_rows; ++i)
      mat(i, i) = static_cast<value_type>(1);
    return mat;
  }

  auto multiRQ() -> matrix_2d
  {
    if (this->m_rows > 3) {
        auto [Q, R] = this->houseHolder_parallel();
        return (R * Q);
    }

    auto [Q, R] = this->houseHolder();
    return (R * Q);
  }

  auto eigenvalues() -> std::vector<value_type>
  {
    std::vector<value_type> eigs(Rows);

    if (Rows == 1) {
        eigs[0] = this->operator()(0, 0);
    } else if (Rows == 2) {
        value_type mean = (this->operator()(0, 0) + this->operator()(1, 1)) / 2;
        value_type d = std::sqrt(mean * mean - this->determinant_laplace_serial());
        eigs[0] = mean - d;
        eigs[1] = mean + d;
    } else {
        matrix_2d B = this->multiRQ();

        for (int n {}; n < 50; ++n) {
            matrix_2d shift = this->identity() * B(this->m_rows - 1, this->m_cols - 1);
            matrix_2d C = B - shift;
            B = C.multiRQ() + shift;
          }

        for (std::size_t i {}; i < this->m_rows; ++i) {
            eigs[i] = (B(i, i));
          }
      }

    return eigs;
  }

  void swap_row(std::size_t row1, std::size_t row2)
  {
    std::vector<value_type> aux(this->m_cols);

    for (std::size_t j {}; j < this->m_cols; ++j) {
        aux[j] = this->operator()(row1, j);
      }
    for (std::size_t j {}; j < this->m_cols; ++j) {
        this->operator()(row1, j) = this->operator()(row2, j);
        this->operator()(row2, j) = aux[j];
      }
  }
  void vec_to_row(std::size_t row, const std::vector<value_type>& vec)
  {
    for (std::size_t j {}; j < this->m_cols; ++j) {
        this->operator()(row, j) = vec[j];
      }
  }

  void vec_to_column(std::size_t col, const std::vector<value_type>& vec)
  {
    for (std::size_t i {}; i < this->m_rows; ++i) {
        this->operator()(i, col) = vec[i];
      }
  }

  auto transpose() -> matrix_2d
  {
    const std::size_t rows = this->m_rows;
    const std::size_t cols = this->m_cols;

    matrix_2d<value_type, rows, cols, ContainerType, AllocatorType> matTrasp;

    if (this->m_rows > std::size_t { 3 } && this->m_cols > std::size_t { 3 }) {
        parallel::parallel_for(parallel::blocked_range { std::size_t {}, this->m_rows },
            [&](auto range)
            {
              for (auto i = range.begin(); i < range.end(); ++i) {
                  for (std::size_t j {}; j < this->m_cols; ++j) {
                      matTrasp(j, i) = this->operator()(i, j);
                    }
                }
            });
    } else {
        for (std::size_t i {}; i < this->m_rows; ++i) {
            for (std::size_t j {}; j < this->m_cols; ++j) {
                matTrasp(j, i) = this->operator()(i, j);
              }
          }
      }

    return matTrasp;
  }

  /// eigenvectors using Gauss Elimination
  auto eigenvectors() -> matrix_2d
  {
    auto eigenvalues = this->eigenvalues();
    const std::size_t num_eigenvalues = eigenvalues.size();

    matrix_2d A = *this;

    auto solve_system = [&](value_type lambda) -> std::vector<value_type>
    {
      // Al = (A - λI)
      matrix_2d Al = A - this->identity() * lambda;

      //  Gauss
      for (std::size_t i {}; i < Al.rows() - 1; ++i) {
          value_type pivot = Al(i, i);
          if (std::abs(pivot) < 1e-10) {
              for (std::size_t n = 1; n < Al.rows(); ++n) {
                  if (i + n < Al.rows()) {
                      if (std::abs(Al(i + 1, i)) > 1.e-10) {
                          Al.swap_row(i, i + 1);
                          pivot = Al(i, i);
                          break;
                      }
                  } else {
                      break;
                    }
                }

              if (std::abs(pivot) < 1e-10) continue;
          }

          for (std::size_t j = i + 1; j < Al.rows(); ++j) {
              value_type factor = Al(j, i) / pivot;
              for (std::size_t k = i; k < Al.columns(); ++k) {
                  Al(j, k) -= factor * Al(i, k);
                }
            }
        }

      std::vector<value_type> eigenvector(Al.rows(), value_type { 1 });

      // get eigenvector

      if (std::abs(Al(Al.rows() - 1, Al.columns() - 1)) > 1e-10)
        eigenvector[Al.rows() - 1] = value_type { 0 };

      for (std::size_t i = Al.rows() - 2; i < Al.rows(); --i) {
          value_type sum {};
          if (std::abs(Al(i, i)) < 1e-10) {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              eigenvector[i + 1] = -sum / Al(i, i + 1);
          } else {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              if (std::abs(sum) < 1e-10) {
                  eigenvector[i] = value_type { 0 };
              } else {
                  eigenvector[i] = -sum / Al(i, i);
                }
            }
        }

#ifdef _MSVC_LANG  // error C2678: '|' binary ; some problem with views
      value_type euclideanNorm = std::sqrt(std::accumulate(eigenvector.begin(), eigenvector.end(),
          value_type {}, [](value_type sum, value_type v) { return sum + v * v; }));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#else
      value_type euclideanNorm = std::sqrt(std::ranges::fold_left(eigenvector
              | std::views::transform([](value_type v) { return v * v; }),
          value_type {}, std::plus<>()));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#endif
      return eigenvector;
    };

    matrix_2d<ElementType, num_eigenvalues, num_eigenvalues> eigenvectors;

    if (num_eigenvalues > 3) {
        parallel::parallel_for(parallel::blocked_range { std::size_t {}, num_eigenvalues },
            [&](auto& range)
            {
              for (std::size_t j = range.begin(); j < range.end(); ++j) {
                  eigenvectors.vec_to_column(j, solve_system(eigenvalues[j]));
                }
            });
    } else {
        for (std::size_t j {}; j < eigenvectors.columns(); ++j) {
            eigenvectors.vec_to_column(j, solve_system(eigenvalues[j]));
          }
      }

    return eigenvectors;
  }

  auto eigenvectors(const std::vector<value_type>& eigenvalues)
  {
    const std::size_t num_eigenvalues = eigenvalues.size();

    matrix_2d A = *this;

    auto solve_system = [&](value_type lambda) -> std::vector<value_type>
    {
      // Al = (A - λI)
      matrix_2d Al = A - this->identity() * lambda;

      //  Gauss
      for (std::size_t i {}; i < Al.rows() - 1; ++i) {
          value_type pivot = Al(i, i);
          if (std::abs(pivot) < 1e-10) {
              for (std::size_t n = 1; n < Al.rows(); ++n) {
                  if (i + n < Al.rows()) {
                      if (std::abs(Al(i + n, i)) > 1.e-10) {
                          Al.swap_row(i, i + n);
                          pivot = Al(i, i);
                          break;
                      }
                  } else {
                      break;
                    }
                }

              if (std::abs(pivot) < 1e-10) continue;
          }

          for (std::size_t j = i + 1; j < Al.rows(); ++j) {
              value_type factor = Al(j, i) / pivot;
              for (std::size_t k = i; k < Al.columns(); ++k) {
                  Al(j, k) -= factor * Al(i, k);
                }
            }
          // std::println("{}", Al);
        }

      std::vector<value_type> eigenvector(Al.rows(), value_type { 1 });

      // get eigenvector

      if (std::abs(Al(Al.rows() - 1, Al.columns() - 1)) > 1e-10)
        eigenvector[Al.rows() - 1] = value_type { 0 };

      for (std::size_t i = Al.rows() - 2; i < Al.rows(); --i) {
          value_type sum {};
          if (std::abs(Al(i, i)) < 1e-10) {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              eigenvector[i + 1] = -sum / Al(i, i + 1);
          } else {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              if (std::abs(sum) < 1e-10) {
                  eigenvector[i] = value_type { 0 };
              } else {
                  eigenvector[i] = -sum / Al(i, i);
                }
            }
          // std::println("lambda = {}, i = {}, eigenvector[i] = {}
          // eigenvector[i+1] = {}", lambda, i, eigenvector[i], eigenvector[i+1]);
        }

#ifdef _MSVC_LANG  // error C2678: '|' binary ; some problem with views
      value_type euclideanNorm = std::sqrt(std::accumulate(eigenvector.begin(), eigenvector.end(),
          value_type {}, [](value_type sum, value_type v) { return sum + v * v; }));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#else
      value_type euclideanNorm = std::sqrt(std::ranges::fold_left(eigenvector
              | std::views::transform([](value_type v) { return v * v; }),
          value_type {}, std::plus<>()));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#endif
      return eigenvector;
    };

    matrix_2d<ElementType, Rows, Cols, ContainerType, AllocatorType> eigenvectors;

    if (num_eigenvalues > 3) {
        parallel::parallel_for(parallel::blocked_range { std::size_t {}, num_eigenvalues },
            [&](auto& range)
            {
              for (std::size_t j = range.begin(); j < range.end(); ++j) {
                  eigenvectors.vec_to_column(j, solve_system(eigenvalues[j]));
                }
            });
    } else {
        for (std::size_t j {}; j < eigenvectors.columns(); ++j) {
            eigenvectors.vec_to_column(j, solve_system(eigenvalues[j]));
          }
      }

    return eigenvectors;
  }

  auto eigenvectors(const value_type& eigenvalue)
  {
    matrix_2d A = *this;

    auto solve_system = [&](value_type lambda) -> std::vector<value_type>
    {
      // Al = (A - λI)
      matrix_2d Al = A - this->identity() * lambda;

      //  Gauss
      for (std::size_t i {}; i < Al.rows() - 1; ++i) {
          value_type pivot = Al(i, i);
          if (std::abs(pivot) < 1e-10) {
              for (std::size_t n = 1; n < Al.rows(); ++n) {
                  if (i + n < Al.rows()) {
                      if (std::abs(Al(i + 1, i)) > 1.e-10) {
                          Al.swap_row(i, i + 1);
                          pivot = Al(i, i);
                          break;
                      }
                  } else {
                      break;
                    }
                }

              if (std::abs(pivot) < 1e-10) continue;
          }

          for (std::size_t j = i + 1; j < Al.rows(); ++j) {
              value_type factor = Al(j, i) / pivot;
              for (std::size_t k = i; k < Al.columns(); ++k) {
                  Al(j, k) -= factor * Al(i, k);
                }
            }
        }

      std::vector<value_type> eigenvector(Al.rows(), value_type { 1 });

      // get eigenvector

      if (std::abs(Al(Al.rows() - 1, Al.columns() - 1)) > 1e-10)
        eigenvector[Al.rows() - 1] = value_type { 0 };

      for (std::size_t i = Al.rows() - 2; i < Al.rows(); --i) {
          value_type sum {};
          if (std::abs(Al(i, i)) < 1e-10) {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              eigenvector[i + 1] = -sum / Al(i, i + 1);
          } else {
              for (std::size_t j = i + 1; j < Al.columns(); ++j) {
                  sum += Al(i, j) * eigenvector[j];
                }
              if (std::abs(sum) < 1e-10) {
                  eigenvector[i] = value_type { 0 };
              } else {
                  eigenvector[i] = -sum / Al(i, i);
                }
            }
        }

#ifdef _MSVC_LANG  // error C2678: '|' binary ; some problem with views
      value_type euclideanNorm = std::sqrt(std::accumulate(eigenvector.begin(), eigenvector.end(),
          value_type {}, [](value_type sum, value_type v) { return sum + v * v; }));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#else
      value_type euclideanNorm = std::sqrt(std::ranges::fold_left(eigenvector
              | std::views::transform([](value_type v) { return v * v; }),
          value_type {}, std::plus<>()));
      std::ranges::for_each(eigenvector, [euclideanNorm](value_type& v) { v /= euclideanNorm; });
#endif
      return eigenvector;
    };

    const std::size_t rows = Rows;
    const std::size_t cols { 1 };
    matrix_2d<ElementType, rows, cols, ContainerType, AllocatorType> eigenvector;
    eigenvector.vec_to_column(std::size_t {}, solve_system(eigenvalue));

    return eigenvector;
  }
};  // class matrix_2d

///   ----  OPERATORS -----

export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
  requires std::is_arithmetic_v<S>
matrix_2d<Element, R, K, Cont, Alloc> operator*(const matrix_2d<Element, R, K, Cont, Alloc>& m,
    S scalar)
{
  auto mm = m;
  mm *= scalar;
  return mm;
}
export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
  requires std::is_arithmetic_v<S>
matrix_2d<Element, R, K, Cont, Alloc> operator*(S scalar,
    const matrix_2d<Element, R, K, Cont, Alloc>& m)
{
  auto mm = m;
  mm *= scalar;
  return mm;
}
export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
  requires std::is_arithmetic_v<S>
matrix_2d<Element, R, K, Cont, Alloc>&& operator*(matrix_2d<Element, R, K, Cont, Alloc>&& m,
    S scalar)
{
  m *= scalar;
  return std::move(m);
}
export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
  requires std::is_arithmetic_v<S>
matrix_2d<Element, R, K, Cont, Alloc>&& operator*(S scalar,
    matrix_2d<Element, R, K, Cont, Alloc>&& m)
{
  m *= scalar;
  return std::move(m);
}
export template<typename Element, std::size_t R, std::size_t K, std::size_t S,
    template<typename...> class Cont, template<typename> class Alloc>
matrix_2d<Element, R, S, Cont, Alloc> operator*(const matrix_2d<Element, R, K, Cont, Alloc>& m1,
    const matrix_2d<Element, K, S, Cont, Alloc>& m2)
{
  if (!(m1.columns() == m2.rows())) {
      throw std::invalid_argument("Cannot multiply matrices.");
  }
  /*
          | a v c|        |q|      |p|
          | d e f|    *   |w|    = |l|
          | g h i|3x3     |j|3x1   |m|3x1
  */
  matrix_2d<Element, R, S, Cont, Alloc> mm;

  parallel::parallel_for(std::size_t {}, m1.rows(),
      [&](auto i)
      {
        for (std::size_t j {}; j < mm.columns(); ++j) {
            for (std::size_t k {}; k < mm.rows(); ++k) {
                mm(i, j) += m1(i, k) * m2(k, j);
              }
          }
      });

  return mm;
}
export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
matrix_2d<Element, R, K, Cont, Alloc> operator/(const matrix_2d<Element, R, K, Cont, Alloc>& m,
    S scalar)
  requires std::is_arithmetic_v<S> && (scalar != 0)
{
  // if (scalar == 0) {
  //     throw std::invalid_argument("Cannot divide by 0.");
  // }

  auto mm = m;
  mm /= scalar;
  return mm;
}
export template<typename Element, std::size_t R, std::size_t K, typename S,
    template<typename...> class Cont, template<typename> class Alloc>
matrix_2d<Element, R, K, Cont, Alloc>&& operator/(matrix_2d<Element, R, K, Cont, Alloc>&& m,
    S scalar)
  requires std::is_arithmetic_v<S> && (scalar != 0)
{
  // if (scalar == 0) {
  //     throw std::invalid_argument("Cannot divide by 0.");
  // }

  m /= scalar;
  return std::move(m);
}
export template<typename Element, std::size_t R, std::size_t K, template<typename...> class Cont,
    template<typename> class Alloc>
matrix_2d<Element, R, K, Cont, Alloc> operator+(const matrix_2d<Element, R, K, Cont, Alloc>& m1,
    const matrix_2d<Element, R, K, Cont, Alloc>& m2)
{
  matrix_2d<Element, R, K, Cont, Alloc> mm;

  parallel::parallel_for(std::size_t {}, m1.rows(),
      [&](auto i)
      {
        for (std::size_t j {}; j < mm.columns(); ++j) {
            mm(i, j) = m1(i, j) + m2(i, j);
          }
      });

  return mm;
}
export template<typename Element, std::size_t R, std::size_t K, template<typename...> class Cont,
    template<typename> class Alloc>
matrix_2d<Element, R, K, Cont, Alloc> operator-(const matrix_2d<Element, R, K, Cont, Alloc>& m1,
    const matrix_2d<Element, R, K, Cont, Alloc>& m2)
{
  matrix_2d<Element, R, K, Cont, Alloc> mm;

  parallel::parallel_for(std::size_t {}, m1.rows(),
      [&](auto i)
      {
        for (std::size_t j {}; j < mm.columns(); ++j) {
            mm(i, j) = m1(i, j) - m2(i, j);
          }
      });

  return mm;
}

#ifdef __clang__
#  ifdef USING_TBBLIB

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using aligned_safe_matrix_2d_t = matrix_2d<ElementType, Rows, Cols, oneapi::tbb::concurrent_vector,
    oneapi::tbb::cache_aligned_allocator>;

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using aligned_fast_matrix_2d_t = matrix_2d<ElementType, Rows, Cols, std::vector,
    oneapi::tbb::cache_aligned_allocator>;

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using scalable_safe_matrix_2d_t = matrix_2d<ElementType, Rows, Cols, oneapi::tbb::concurrent_vector,
    oneapi::tbb::scalable_allocator>;

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using scalable_fast_matrix_2d_t = matrix_2d<ElementType, Rows, Cols, std::vector,
    oneapi::tbb::scalable_allocator>;

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using mat = matrix_2d<ElementType, Rows, Cols, std::vector, oneapi::tbb::scalable_allocator>;

export template<std::size_t Rows, std::size_t Cols>
using fmat = matrix_2d<float, Rows, Cols, std::vector, oneapi::tbb::scalable_allocator>;

export template<std::size_t Rows, std::size_t Cols>
using dmat = matrix_2d<double, Rows, Cols, std::vector, oneapi::tbb::scalable_allocator>;

export template<std::size_t Rows, std::size_t Cols>
using imat = matrix_2d<int, Rows, Cols, std::vector, oneapi::tbb::scalable_allocator>;
#  else

export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using mat = matrix_2d<ElementType, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using fmat = matrix_2d<float, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using dmat = matrix_2d<double, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using imat = matrix_2d<int, Rows, Cols, std::vector>;

#  endif
#else
export template<typename ElementType, std::size_t Rows, std::size_t Cols>
using mat = matrix_2d<ElementType, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using fmat = matrix_2d<float, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using dmat = matrix_2d<double, Rows, Cols, std::vector>;

export template<std::size_t Rows, std::size_t Cols>
using imat = matrix_2d<int, Rows, Cols, std::vector>;

#endif
}  // namespace jf::matrix

#ifdef USING_FMTLIB
template<typename ElementType, std::size_t Rows, std::size_t Cols,
    template<typename...> class ContainerType, template<typename> class AllocatorType>
struct fmt::formatter<jf::matrix::matrix_2d<ElementType, Rows, Cols, ContainerType, AllocatorType>>
{
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const jf::matrix::matrix_2d<ElementType, Rows, Cols, ContainerType, AllocatorType>& m,
      FormatContext& ctx) const
  {
    auto out = ctx.out();
    if (m.empty()) {
        return fmt::format_to(out, "[ ]");
    }

    fmt::format_to(out, "\n");
    for (std::size_t i {}; i < m.rows(); ++i) {
        fmt::format_to(out, "[");
        for (std::size_t j {}; j < m.columns(); ++j) {
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
#else

template<typename ElementType, std::size_t Rows, std::size_t Cols,
    template<typename...> class ContainerType, template<typename> class AllocatorType>
struct std::formatter<jf::matrix::matrix_2d<ElementType, Rows, Cols, ContainerType, AllocatorType>>
{
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const jf::matrix::matrix_2d<ElementType, Rows, Cols, ContainerType, AllocatorType>& m,
      FormatContext& ctx) const
  {
    auto out = ctx.out();
    if (m.empty()) {
        return std::format_to(out, "[ ]");
    }

    std::format_to(out, "\n");
    for (std::size_t i {}; i < m.rows(); ++i) {
        std::format_to(out, "[");
        for (std::size_t j {}; j < m.columns(); ++j) {
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
#endif
