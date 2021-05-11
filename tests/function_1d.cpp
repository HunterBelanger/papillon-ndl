#include <gtest/gtest.h>

#include <PapillonNDL/difference_1d.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/sum_1d.hpp>
#include <memory>
#include <vector>

namespace pndl {
namespace {

//==============================================================================
// Region1D
TEST(Region1D, MinMaxX) {
  std::vector<double> x{1., 2., 3., 4., 5., 6.};
  std::vector<double> y{1., 2., 3., 4., 5., 6.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x, y, interp);

  EXPECT_DOUBLE_EQ(R.min_x(), x.front());
  EXPECT_DOUBLE_EQ(R.max_x(), x.back());
}

TEST(Region1D, Interpolation) {
  std::vector<double> x{1., 2., 3., 4., 5., 6.};
  std::vector<double> y{1., 2., 3., 4., 5., 6.};
  Interpolation interp = Interpolation::LogLin;

  Region1D R(x, y, interp);

  EXPECT_EQ(R.interpolation()[0], interp);
}

TEST(Region1D, EvaluationContinuous) {
  std::vector<double> x_vals{1., 2., 3., 4., 5., 6.};
  std::vector<double> y_vals{1., 2., 3., 4., 5., 6.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x_vals, y_vals, interp);

  EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front() - 0.001));

  EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()));

  double x = 1.23456;
  EXPECT_DOUBLE_EQ(x, R(x));

  x = 3.238197263;
  EXPECT_DOUBLE_EQ(x, R(x));

  x = 4.238197263;
  EXPECT_DOUBLE_EQ(x, R(x));

  EXPECT_DOUBLE_EQ(x_vals.back(), R(x_vals.back()));

  EXPECT_DOUBLE_EQ(x_vals.back(), R(x_vals.back() + 0.001));
}

TEST(Region1D, EvaluationDiscontinuous) {
  std::vector<double> x_vals{1., 2., 3., 3., 4., 5.};
  std::vector<double> y_vals{1., 2., 3., 6., 8., 10.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x_vals, y_vals, interp);

  EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front() - 0.001));

  EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()));

  double x = 1.23456;
  EXPECT_DOUBLE_EQ(x, R(x));

  x = 2.999999999;
  EXPECT_DOUBLE_EQ(x, R(x));

  x = 3.000000001;
  EXPECT_DOUBLE_EQ(2 * x, R(x));

  x = 4.238197263;
  EXPECT_DOUBLE_EQ(2 * x, R(x));

  EXPECT_DOUBLE_EQ(2 * x_vals.back(), R(x_vals.back()));

  EXPECT_DOUBLE_EQ(2 * x_vals.back(), R(x_vals.back() + 0.001));
}

TEST(Region1D, Integration) {
  std::vector<double> x_vals{0., 2., 4., 6.};
  std::vector<double> y_vals{1., 1., 3., 2.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x_vals, y_vals, interp);

  double x_low = 0.;
  double x_hi = 2.;
  double i1 = 2.;
  EXPECT_DOUBLE_EQ(i1, R.integrate(x_low, x_hi));

  x_low = 2.;
  x_hi = 4.;
  double i2 = 4.;
  EXPECT_DOUBLE_EQ(i2, R.integrate(x_low, x_hi));

  x_low = 4.;
  x_hi = 6.;
  double i3 = 5.;
  EXPECT_DOUBLE_EQ(i3, R.integrate(x_low, x_hi));

  x_low = 0.;
  x_hi = 6.;
  EXPECT_DOUBLE_EQ(i1 + i2 + i3, R.integrate(x_low, x_hi));

  x_low = 3.;
  x_hi = 5.;
  double i = 5.25;
  EXPECT_DOUBLE_EQ(i, R.integrate(x_low, x_hi));

  // Ensure integration with inverted limits comes out as the negative
  EXPECT_DOUBLE_EQ(-i, R.integrate(x_hi, x_low));
}

TEST(Region1D, Size) {
  std::vector<double> x_vals{1., 2., 3., 4., 5., 6.};
  std::vector<double> y_vals{1., 2., 3., 4., 5., 6.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x_vals, y_vals, interp);

  EXPECT_EQ(x_vals.size(), R.size());
}

TEST(Region1D, XY) {
  std::vector<double> x_vals{1., 2., 3., 4., 5., 6.};
  std::vector<double> y_vals{1., 2., 3., 4., 5., 6.};
  Interpolation interp = Interpolation::LinLin;

  Region1D R(x_vals, y_vals, interp);

  const std::vector<double>& xref = R.x();
  const std::vector<double>& yref = R.y();

  EXPECT_EQ(xref.size(), yref.size());
  EXPECT_EQ(xref.size(), R.size());
  EXPECT_EQ(x_vals.size(), y_vals.size());
  EXPECT_EQ(x_vals.size(), xref.size());

  for (size_t i = 0; i < xref.size(); i++) {
    EXPECT_DOUBLE_EQ(xref[i], x_vals[i]);
    EXPECT_DOUBLE_EQ(yref[i], y_vals[i]);
  }
}

//==============================================================================
// MultiRegion1D
TEST(MultiRegion1D, ConstructorRegions) {
  std::vector<Region1D> regions{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({5., 6.}, {2., 23.}, Interpolation::LinLin)};
  EXPECT_ANY_THROW(MultiRegion1D m_fail(regions));

  std::vector<Region1D> regions2{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({-1., 6.}, {2., 23.}, Interpolation::LinLin)};
  EXPECT_ANY_THROW(MultiRegion1D m_fail2(regions2));

  std::vector<Region1D> regions3{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({4., 6.}, {2., 23.}, Interpolation::LinLin)};
  EXPECT_NO_THROW(MultiRegion1D m(regions3));
}

TEST(MultiRegion1D, ConstructorACE) {
  std::vector<uint32_t> NBT1{2, 4, 5};
  std::vector<Interpolation> INT1{Interpolation::LinLin, Interpolation::LinLin,
                                  Interpolation::LinLin};
  std::vector<double> x1{1., 2., 6., 10.};
  std::vector<double> y1{0., 1., 2., 6., 20.};
  EXPECT_ANY_THROW(MultiRegion1D R(NBT1, INT1, x1, y1));

  std::vector<uint32_t> NBT2{2, 5};
  std::vector<Interpolation> INT2{Interpolation::LinLin, Interpolation::LinLin,
                                  Interpolation::LinLin};
  std::vector<double> x2{1., 2., 2., 6., 10.};
  std::vector<double> y2{0., 1., 2., 6., 20.};
  EXPECT_ANY_THROW(MultiRegion1D R(NBT2, INT2, x2, y2));

  std::vector<uint32_t> NBT3{2, 4, 5};
  std::vector<Interpolation> INT3{Interpolation::LinLin, Interpolation::LinLin,
                                  Interpolation::LinLin};
  std::vector<double> x3{1., 2., 2., 1.5, 10.};
  std::vector<double> y3{0., 1., 2., 6., 20.};
  EXPECT_ANY_THROW(MultiRegion1D R(NBT3, INT3, x3, y3));

  std::vector<uint32_t> NBT{2, 4, 5};
  std::vector<Interpolation> INT{Interpolation::LinLin, Interpolation::LinLin,
                                 Interpolation::LinLin};
  std::vector<double> x{1., 2., 2., 6., 10.};
  std::vector<double> y{0., 1., 2., 6., 20.};
  EXPECT_NO_THROW(MultiRegion1D R(NBT, INT, x, y));
}

TEST(MultiRegion1D, Size) {
  std::vector<Region1D> regions{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({4., 6.}, {2., 23.}, Interpolation::LinLin)};
  MultiRegion1D R(regions);

  EXPECT_EQ(R.size(), regions.size());
}

TEST(MultiRegion1D, Region) {
  std::vector<Region1D> regions{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({4., 6.}, {2., 23.}, Interpolation::LinLin)};
  MultiRegion1D R(regions);

  EXPECT_EQ(R.size(), regions.size());
}

TEST(MultiRegion1D, MinMaxX) {
  std::vector<Region1D> regions{
      Region1D({1., 4.}, {1., 2.}, Interpolation::LinLin),
      Region1D({4., 6.}, {2., 23.}, Interpolation::LinLin)};
  MultiRegion1D R(regions);

  EXPECT_DOUBLE_EQ(regions.front().min_x(), R.min_x());
  EXPECT_DOUBLE_EQ(regions.back().max_x(), R.max_x());
}

TEST(MultiRegion1D, Evaluation) {
  std::vector<Region1D> regions{
      Region1D({1., 4.}, {1., 4.}, Interpolation::LinLin),
      Region1D({4., 8.}, {8., 16.}, Interpolation::LinLin)};
  MultiRegion1D R(regions);

  double x = 0.9;
  double y = 1.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 1.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 3.;
  y = 3.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 3.999999;
  y = 3.999999;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 4.0000001;
  y = 2 * x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 7.;
  y = 2 * x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 8.;
  y = 2 * x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 12.;
  EXPECT_DOUBLE_EQ(y, R(x));
}

TEST(MultiRegion1D, EvaluationDiscontinuous) {
  std::vector<uint32_t> NBT{2, 4, 5};
  std::vector<Interpolation> INT{Interpolation::LinLin, Interpolation::LinLin,
                                 Interpolation::LinLin};
  std::vector<double> x_vals{1., 2., 2., 6., 10.};
  std::vector<double> y_vals{0., 1., 2., 6., 20.};
  MultiRegion1D R(NBT, INT, x_vals, y_vals);

  double x = x_vals.front();
  double y = y_vals.front();
  EXPECT_DOUBLE_EQ(y, R(x));

  x -= 1.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 1.99999999;
  y = x - 1.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 2.000001;
  y = x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 3.;
  y = x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 5.9999999999;
  y = x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 6.;
  y = x;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 9.;
  y = 16.5;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 10.;
  y = 20.;
  EXPECT_DOUBLE_EQ(y, R(x));

  x = 15.;
  EXPECT_DOUBLE_EQ(y, R(x));
}

TEST(MultiRegion1D, Integration) {
  std::vector<uint32_t> NBT{2, 3, 4};
  std::vector<Interpolation> INT{Interpolation::LinLin, Interpolation::LinLin,
                                 Interpolation::LinLin};
  std::vector<double> x_vals{0., 2., 4., 6.};
  std::vector<double> y_vals{1., 1., 3., 2.};

  MultiRegion1D R(NBT, INT, x_vals, y_vals);
  EXPECT_EQ(3, R.size());

  double x_low = 0.;
  double x_hi = 2.;
  double i1 = 2.;
  EXPECT_DOUBLE_EQ(i1, R.integrate(x_low, x_hi));

  x_low = 2.;
  x_hi = 4.;
  double i2 = 4.;
  EXPECT_DOUBLE_EQ(i2, R.integrate(x_low, x_hi));

  x_low = 4.;
  x_hi = 6.;
  double i3 = 5.;
  EXPECT_DOUBLE_EQ(i3, R.integrate(x_low, x_hi));

  x_low = 0.;
  x_hi = 6.;
  EXPECT_DOUBLE_EQ(i1 + i2 + i3, R.integrate(x_low, x_hi));

  x_low = 3.;
  x_hi = 5.;
  double i = 5.25;
  EXPECT_DOUBLE_EQ(i, R.integrate(x_low, x_hi));

  // Ensure integration with inverted limits produces negative
  EXPECT_DOUBLE_EQ(-i, R.integrate(x_hi, x_low));
}

//==============================================================================
// Polynomial1D
TEST(Polynomial1D, Order) {
  std::vector<double> coeffs{3., 4., 5., 6.};
  Polynomial1D poly(coeffs);

  EXPECT_EQ(poly.order(), coeffs.size() - 1);

  std::vector<double> coeffs2{3., 4., 5., 6., 2., 1., 3.5, 6.5};
  Polynomial1D poly2(coeffs2);

  EXPECT_EQ(poly2.order(), coeffs2.size() - 1);
}

TEST(Polynomial1D, Coefficients) {
  std::vector<double> coeffs{1.1, 2.2, 3.3, 4.4};
  Polynomial1D poly(coeffs);

  for (size_t i = 0; i <= poly.order(); i++) {
    EXPECT_DOUBLE_EQ(coeffs[i], poly.coefficient(i));
  }
}

TEST(Polynomial1D, Evaluation) {
  std::vector<double> coeffs{1.1, 2.2, 3.3, 4.4};
  Polynomial1D poly(coeffs);

  double x = 0.;
  double y = coeffs[0];
  EXPECT_DOUBLE_EQ(y, poly(x));

  x = 1.;
  y = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3];
  EXPECT_DOUBLE_EQ(y, poly(x));

  x = 2.;
  y = 53.900000000000006;
  EXPECT_DOUBLE_EQ(y, poly(x));

  x = 5.;
  y = 644.6;
  EXPECT_DOUBLE_EQ(y, poly(x));

  x = 20.;
  y = 36565.1;
  EXPECT_DOUBLE_EQ(y, poly(x));
}

TEST(Polynomial1D, Integration) {
  std::vector<double> coeffs{1.1, 2.2, 3.3, 4.4};
  Polynomial1D poly(coeffs);

  double x_low = 1.;
  double x_hi = 5.;
  double i = 853.6;
  EXPECT_DOUBLE_EQ(i, poly.integrate(x_low, x_hi));

  x_low = -7.8;
  x_hi = 22.7;
  i = 301926.74985;
  EXPECT_DOUBLE_EQ(i, poly.integrate(x_low, x_hi));

  x_hi = -7.8;
  x_low = 22.7;
  EXPECT_DOUBLE_EQ(-i, poly.integrate(x_low, x_hi));
}

//==============================================================================
// Sum1D
TEST(Sum1D, Evaluation) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Sum1D sum(t1, t2);

  EXPECT_DOUBLE_EQ(sum(1.), t1->evaluate(1.) + t2->evaluate(1.));
  EXPECT_DOUBLE_EQ(sum(1.5), t1->evaluate(1.5) + t2->evaluate(1.5));
  EXPECT_DOUBLE_EQ(sum(3.), t1->evaluate(3.) + t2->evaluate(3.));
}

TEST(Sum1D, Integration) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Sum1D sum(t1, t2);

  EXPECT_DOUBLE_EQ(sum.integrate(1., 3.),
                   t1->integrate(1., 3.) + t2->integrate(1., 3.));
  EXPECT_DOUBLE_EQ(sum.integrate(1.5, 2.5),
                   t1->integrate(1.5, 2.5) + t2->integrate(1.5, 2.5));
  EXPECT_DOUBLE_EQ(sum.integrate(2.5, 1.5),
                   t1->integrate(2.5, 1.5) + t2->integrate(2.5, 1.5));
}

TEST(Sum1D, Terms) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Sum1D sum(t1, t2);

  EXPECT_EQ(sum.term_1().shared_from_this(), t1);
  EXPECT_EQ(sum.term_2().shared_from_this(), t2);
}

//==============================================================================
// Difference1D
TEST(Difference1D, Evaluation) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Difference1D diff(t1, t2);

  EXPECT_DOUBLE_EQ(diff(1.), t1->evaluate(1.) - t2->evaluate(1.));
  EXPECT_DOUBLE_EQ(diff(1.5), t1->evaluate(1.5) - t2->evaluate(1.5));
  EXPECT_DOUBLE_EQ(diff(3.), t1->evaluate(3.) - t2->evaluate(3.));
}

TEST(Difference1D, Integration) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Difference1D diff(t1, t2);

  EXPECT_DOUBLE_EQ(diff.integrate(1., 3.),
                   t1->integrate(1., 3.) - t2->integrate(1., 3.));
  EXPECT_DOUBLE_EQ(diff.integrate(1.5, 2.5),
                   t1->integrate(1.5, 2.5) - t2->integrate(1.5, 2.5));
  EXPECT_DOUBLE_EQ(diff.integrate(2.5, 1.5),
                   t1->integrate(2.5, 1.5) - t2->integrate(2.5, 1.5));
}

TEST(Difference1D, Terms) {
  std::shared_ptr<Region1D> t1 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 4., 9.}),
      Interpolation::LinLin);
  std::shared_ptr<Region1D> t2 = std::make_shared<Region1D>(
      std::vector<double>({1., 2., 3.}), std::vector<double>({1., 2., 3.}),
      Interpolation::LinLin);
  Difference1D diff(t1, t2);

  EXPECT_EQ(diff.term_1().shared_from_this(), t1);
  EXPECT_EQ(diff.term_2().shared_from_this(), t2);
}

}  // namespace
}  // namespace pndl