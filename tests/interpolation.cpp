#include <gtest/gtest.h>

#include <PapillonNDL/interpolation.hpp>
#include <vector>

namespace pndl {
namespace {

//==============================================================================
// Histogram Tests
TEST(Histogram, Interpolate) {
  Histogram interp;

  double x1 = 8.;
  double y1 = 4.;

  double x2 = 10.;
  double y2 = 6.;

  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 9.;
  double y = y1;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  x = x2 - 0.000000001;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(Histogram, Invert) {
  Histogram interp;

  double x1 = 8.;
  double y1 = 4.;

  double x2 = 10.;
  double y2 = 6.;

  EXPECT_DOUBLE_EQ(x1, interp.invert(x1, x1, y1, x2, y2));

  double x = 9.;
  double y = y1;
  EXPECT_DOUBLE_EQ(x1, interp.invert(x, x1, y1, x2, y2));

  x = x2 - 0.000000001;
  EXPECT_DOUBLE_EQ(x1, interp.invert(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x1, interp.invert(x2, x1, y1, x2, y2));
}

TEST(Histogram, Integrate) {
  Histogram interp;
  double x1 = 0.;
  double y1 = 4.;
  double x2 = 5.;
  double y2 = 8.;

  double x_low = x1;
  double x_hi = x2;
  double i = y1 * (x_hi - x_low);
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 1.3;
  x_hi = 4.7;
  i = y1 * (x_hi - x_low);
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(Histogram, VerifyXGrid) {
  Histogram interp;

  std::vector<double> x_bad{1., 2., 2., 4.};
  std::vector<double> x_good{
      1.,
      2.,
      3.,
      4.,
  };

  EXPECT_THROW(interp.verify_x_grid(x_bad.begin(), x_bad.end()), PNDLException);
  EXPECT_NO_THROW(interp.verify_x_grid(x_good.begin(), x_good.end()));
}

TEST(Histogram, VerifyYGrid) {
  Histogram interp;
  std::vector<double> y{-2., 1., 3., 3., -7.};
  EXPECT_NO_THROW(interp.verify_y_grid(y.begin(), y.end()));
}

//==============================================================================
// LinLin Tests
TEST(LinLin, Interpolate) {
  LinLin interp;
  // Here, y is linear in x
  // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1

  double x1 = 8.;
  double y1 = 0.;

  double x2 = 10.;
  double y2 = 6.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 9.;
  double y = 3.;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(LinLin, Invert) {
  LinLin interp;
  // Here, y is linear in x
  // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1

  double x1 = 8.;
  double y1 = 0.;

  double x2 = 10.;
  double y2 = 6.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  double x = 9.;
  double y = 3.;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));
}

TEST(LinLin, Integrate) {
  LinLin interp;
  double x1 = 0.;
  double y1 = 0.;
  double x2 = 5.;
  double y2 = 5.;

  double x_low = x1;
  double x_hi = x2;
  double i = 12.5;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 1.3;
  x_hi = 4.7;
  i = 10.2;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(LinLin, VerifyXGrid) {
  LinLin interp;

  std::vector<double> x_good{-2., -1., 4., 8.};
  std::vector<double> x_bad{-5., 0., -1., 6.};

  EXPECT_NO_THROW(interp.verify_x_grid(x_good.begin(), x_good.end()));
  EXPECT_THROW(interp.verify_x_grid(x_bad.begin(), x_bad.end()), PNDLException);
}

TEST(LinLin, VerifyYGrid) {
  LinLin interp;

  std::vector<double> y1{-2., -1., 4., 8.};
  std::vector<double> y2{-5., 0., -1., 6.};

  EXPECT_NO_THROW(interp.verify_y_grid(y1.begin(), y1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y2.begin(), y2.end()));
}

//==============================================================================
// LogLin Tests
TEST(LogLin, Interpolate) {
  // Here, y is linear in log(x)
  // log(y) = ((x - x1)/(x2 - x1))*log(y2/y1) + log(y1)
  LogLin interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.4422495703074083;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  x = 3.;
  y = 2.080083823051904;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(LogLin, Invert) {
  LogLin interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.4422495703074083;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  x = 3.;
  y = 2.080083823051904;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));
}

TEST(LogLin, Integrate) {
  LogLin interp;
  double x1 = -10.;
  double y1 = 1.;
  double x2 = 15.;
  double y2 = 11.;

  double x_low = x1;
  double x_hi = x2;
  double i = 104.25809785606157;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = -4.378;
  x_hi = 5.69;
  i = 29.078502985642945;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(LogLin, VerifyXGrid) {
  LogLin interp;
  std::vector<double> x1{-5., -4., -3., -2.};
  std::vector<double> x2{2., 3., 4., 5.};
  std::vector<double> x3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_x_grid(x1.begin(), x1.end()));
  EXPECT_NO_THROW(interp.verify_x_grid(x2.begin(), x2.end()));
  EXPECT_NO_THROW(interp.verify_x_grid(x3.begin(), x3.end()));
}

TEST(LogLin, VerifyYGrid) {
  LogLin interp;
  std::vector<double> y_good_1{-5., -4., -3., -2.};
  std::vector<double> y_good_2{2., 3., 4., 5.};
  std::vector<double> y_bad{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_y_grid(y_good_1.begin(), y_good_1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y_good_2.begin(), y_good_2.end()));
  EXPECT_THROW(interp.verify_y_grid(y_bad.begin(), y_bad.end()), PNDLException);
}

//==============================================================================
// LinLog Tests
TEST(LinLog, Interpolate) {
  // Here, y is linear in log(x)
  // y = (log(x/x1)/log(x2/x1))*(y2 - y1) + y1
  LinLog interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 2.;
  double y = 2.;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  x = 3.;
  y = 2.584962500721156;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(LinLog, Invert) {
  LinLog interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  double x = 2.;
  double y = 2.;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  x = 3.;
  y = 2.584962500721156;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));
}

TEST(LinLog, Integrate) {
  LinLog interp;
  double x1 = 5.;
  double y1 = 1.;
  double x2 = 15.;
  double y2 = 21.;

  double x_low = x1;
  double x_hi = x2;
  double i = 127.95215467463257;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 6.4;
  x_hi = 12.8;
  i = 80.169216888001344;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(LinLog, VerifyXGrid) {
  LinLog interp;
  std::vector<double> x_good_1{-5., -4., -3., -2.};
  std::vector<double> x_good_2{2., 3., 4., 5.};
  std::vector<double> x_bad{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_x_grid(x_good_1.begin(), x_good_1.end()));
  EXPECT_NO_THROW(interp.verify_x_grid(x_good_2.begin(), x_good_2.end()));
  EXPECT_THROW(interp.verify_x_grid(x_bad.begin(), x_bad.end()), PNDLException);
}

TEST(LinLog, VerifyYGrid) {
  LinLog interp;
  std::vector<double> y1{-5., -4., -3., -2.};
  std::vector<double> y2{2., 3., 4., 5.};
  std::vector<double> y3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_y_grid(y1.begin(), y1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y2.begin(), y2.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y3.begin(), y3.end()));
}

//==============================================================================
// LogLog Tests
TEST(LogLog, Interpolate) {
  // Here, log(y) is linear in log(x)
  // log(y) = (log(x/x1)/log(x2/x1))*log(y2/y1) + log(y1)
  LogLog interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.7320508075688774;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  x = 3.;
  y = 2.388414221757005;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(LogLog, Invert) {
  LogLog interp;

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.7320508075688774;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  x = 3.;
  y = 2.388414221757005;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));
}

TEST(LogLog, Integrate) {
  LogLog interp;
  double x1 = 1;
  double y1 = 5.;
  double x2 = 5.;
  double y2 = 8.;

  double x_low = x1;
  double x_hi = x2;
  double i = 27.089161107019226;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 2.348;
  x_hi = 3.892;
  i = 10.739759942209645;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(LogLog, VerifyXGrid) {
  LogLog interp;
  std::vector<double> x1{-5., -4., -3., -2.};
  std::vector<double> x2{2., 3., 4., 5.};
  std::vector<double> x3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_x_grid(x1.begin(), x1.end()));
  EXPECT_NO_THROW(interp.verify_x_grid(x2.begin(), x2.end()));
  EXPECT_THROW(interp.verify_y_grid(x3.begin(), x3.end()), PNDLException);
}

TEST(LogLog, VerifyYGrid) {
  LogLog interp;
  std::vector<double> y1{-5., -4., -3., -2.};
  std::vector<double> y2{2., 3., 4., 5.};
  std::vector<double> y3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_y_grid(y1.begin(), y1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y2.begin(), y2.end()));
  EXPECT_THROW(interp.verify_y_grid(y3.begin(), y3.end()), PNDLException);
}

//==============================================================================
// Interpolator Tests
TEST(Interpolator, Interpolate) {
  // Here, log(y) is linear in log(x)
  // log(y) = (log(x/x1)/log(x2/x1))*log(y2/y1) + log(y1)
  Interpolator interp(Interpolation::LogLog);

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.7320508075688774;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  x = 3.;
  y = 2.388414221757005;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));

  interp = Interpolator(Interpolation::LinLin);
  x1 = 8.;
  y1 = 0.;

  x2 = 10.;
  y2 = 6.;
  EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

  x = 9.;
  y = 3.;
  EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
}

TEST(Interpolator, Invert) {
  Interpolator interp(Interpolation::LogLog);

  double x1 = 1.;
  double y1 = 1.;

  double x2 = 4.;
  double y2 = 3.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  double x = 2.;
  double y = 1.7320508075688774;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  x = 3.;
  y = 2.388414221757005;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));

  interp = Interpolator(Interpolation::LinLin);
  x1 = 8.;
  y1 = 0.;

  x2 = 10.;
  y2 = 6.;
  EXPECT_DOUBLE_EQ(x1, interp.invert(y1, x1, y1, x2, y2));

  x = 9.;
  y = 3.;
  EXPECT_DOUBLE_EQ(x, interp.invert(y, x1, y1, x2, y2));

  EXPECT_DOUBLE_EQ(x2, interp.invert(y2, x1, y1, x2, y2));
}

TEST(Interpolator, Integrate) {
  Interpolator interp(Interpolation::LogLog);
  double x1 = 1;
  double y1 = 5.;
  double x2 = 5.;
  double y2 = 8.;

  double x_low = x1;
  double x_hi = x2;
  double i = 27.089161107019226;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 2.348;
  x_hi = 3.892;
  i = 10.739759942209645;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  interp = Interpolator(Interpolation::LinLin);
  x1 = 0.;
  y1 = 0.;
  x2 = 5.;
  y2 = 5.;

  x_low = x1;
  x_hi = x2;
  i = 12.5;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = x2;
  x_hi = x1;
  EXPECT_DOUBLE_EQ(-i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));

  x_low = 1.3;
  x_hi = 4.7;
  i = 10.2;
  EXPECT_DOUBLE_EQ(i, interp.integrate(x_low, x_hi, x1, y1, x2, y2));
}

TEST(Interpolator, VerifyXGrid) {
  Interpolator interp(Interpolation::LogLog);
  std::vector<double> x1{-5., -4., -3., -2.};
  std::vector<double> x2{2., 3., 4., 5.};
  std::vector<double> x3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_x_grid(x1.begin(), x1.end()));
  EXPECT_NO_THROW(interp.verify_x_grid(x2.begin(), x2.end()));
  EXPECT_THROW(interp.verify_y_grid(x3.begin(), x3.end()), PNDLException);

  interp = Interpolator(Interpolation::LinLin);
  std::vector<double> x_good{-2., -1., 4., 8.};
  std::vector<double> x_bad{-5., 0., -1., 6.};

  EXPECT_NO_THROW(interp.verify_x_grid(x_good.begin(), x_good.end()));
  EXPECT_THROW(interp.verify_x_grid(x_bad.begin(), x_bad.end()), PNDLException);
}

TEST(Interpolator, VerifyYGrid) {
  Interpolator interp(Interpolation::LogLog);
  std::vector<double> y1{-5., -4., -3., -2.};
  std::vector<double> y2{2., 3., 4., 5.};
  std::vector<double> y3{-2., -1., 4., 5.};

  EXPECT_NO_THROW(interp.verify_y_grid(y1.begin(), y1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y2.begin(), y2.end()));
  EXPECT_THROW(interp.verify_y_grid(y3.begin(), y3.end()), PNDLException);

  interp = Interpolator(Interpolation::LinLin);
  y1 = {-2., -1., 4., 8.};
  y2 = {-5., 0., -1., 6.};

  EXPECT_NO_THROW(interp.verify_y_grid(y1.begin(), y1.end()));
  EXPECT_NO_THROW(interp.verify_y_grid(y2.begin(), y2.end()));
}

}  // namespace
}  // namespace pndl
