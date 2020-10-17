/*
 * Copyright 2020, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#ifndef PAPILLON_NDL_SHARED_SPAN_H
#define PAPILLON_NDL_SHARED_SPAN_H

#include <iterator>
#include <memory>
#include <vector>

namespace pndl {

template <class T>
class shared_span {
 public:
  // constants and types
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = element_type*;
  using const_pointer = const element_type*;
  using reference = element_type&;
  using const_reference = const element_type&;
  using iterator = pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;

  // constructors, copy, move, and assignment
  template <class InputIt>
  shared_span(InputIt first, InputIt last);
  shared_span(std::initializer_list<element_type> init);
  shared_span(const shared_span& other, std::size_t Offset, std::size_t Count);

  shared_span(const shared_span& other) = default;
  shared_span(shared_span&& other) = default;
  constexpr shared_span& operator=(const shared_span& other) noexcept = default;
  ~shared_span() noexcept = default;

  // subviews
  constexpr shared_span first(std::size_t Count) const;
  constexpr shared_span last(std::size_t Count) const;
  constexpr shared_span subspan(std::size_t Offset, std::size_t Count) const;

  // observers
  constexpr size_type size() const noexcept;
  constexpr size_type size_bytes() const noexcept;
  constexpr bool empty() const noexcept;
  long use_count() const noexcept;

  // element access
  constexpr reference operator[](size_type idx) const;
  constexpr reference front() const;
  constexpr reference back() const;
  constexpr pointer data() const noexcept;

  // iterator support
  constexpr iterator begin() const noexcept;
  constexpr iterator end() const noexcept;
  constexpr reverse_iterator rbegin() const noexcept;
  constexpr reverse_iterator rend() const noexcept;

 private:
  std::shared_ptr<std::vector<element_type>> data_;
  size_type begin_;
  size_type end_;
};

template <class T>
inline shared_span<T>::shared_span(std::initializer_list<element_type> init)
    : data_{nullptr}, begin_{0}, end_{0} {
  data_ = std::make_shared<std::vector<element_type>>(init);
  begin_ = 0;
  end_ = data_->size();
}

template <class T>
template <class InputIt>
inline shared_span<T>::shared_span(InputIt first, InputIt last)
    : data_{nullptr}, begin_{0}, end_{0} {
  data_ = std::make_shared<std::vector<element_type>>(first, last);
  begin_ = 0;
  end_ = data_->size();
}

template <class T>
inline shared_span<T>::shared_span(const shared_span& other, std::size_t Offset,
                                   std::size_t Count)
    : data_{other.data_}, begin_{other.begin_}, end_{other.end_} {
  begin_ += Offset;
  if (begin_ + Count <= end_) end_ = begin_ + Count;
}

template <class T>
inline constexpr shared_span<T> shared_span<T>::first(std::size_t Count) const {
  return shared_span(*this, 0, Count);
}

template <class T>
inline constexpr shared_span<T> shared_span<T>::last(std::size_t Count) const {
  if (Count > (end_ - begin_)) return shared_span(*this);
  return shared_span(*this, end_ - Count, end_);
}

template <class T>
inline constexpr shared_span<T> shared_span<T>::subspan(
    std::size_t Offset, std::size_t Count) const {
  return {*this, Offset, Count};
}

template <class T>
inline constexpr typename shared_span<T>::size_type shared_span<T>::size()
    const noexcept {
  return end_ - begin_;
}

template <class T>
inline constexpr typename shared_span<T>::size_type shared_span<T>::size_bytes()
    const noexcept {
  return (end_ - begin_) * sizeof(element_type);
}

template <class T>
inline constexpr bool shared_span<T>::empty() const noexcept {
  return (end_ - begin_) == 0;
}

template <class T>
inline long shared_span<T>::use_count() const noexcept {
  return data_.use_count();
}

template <class T>
inline constexpr typename shared_span<T>::reference shared_span<T>::operator[](
    size_type idx) const {
  return (*data_)[begin_ + idx];
}

template <class T>
inline constexpr typename shared_span<T>::reference shared_span<T>::front()
    const {
  return (*data_)[begin_];
}

template <class T>
inline constexpr typename shared_span<T>::reference shared_span<T>::back()
    const {
  return (*data_)[end_ - 1];
}

template <class T>
inline constexpr typename shared_span<T>::pointer shared_span<T>::data()
    const noexcept {
  return data_->data() + begin_;
}

template <class T>
inline constexpr typename shared_span<T>::iterator shared_span<T>::begin()
    const noexcept {
  return data_->data() + begin_;
}

template <class T>
inline constexpr typename shared_span<T>::iterator shared_span<T>::end()
    const noexcept {
  return data_->data() + end_;
}

template <class T>
inline constexpr typename shared_span<T>::reverse_iterator
shared_span<T>::rbegin() const noexcept {
  return reverse_iterator(data_->data() + end_);
}

template <class T>
inline constexpr typename shared_span<T>::reverse_iterator
shared_span<T>::rend() const noexcept {
  return reverse_iterator(data_->data() + begin_);
}

}  // namespace pndl

#endif
