#include <PapillonNDL/integration.hpp>
#include <gtest/gtest.h>

namespace {
  using namespace pndl;

  TEST(Integration, Histogram) {
    Interpolation interp = Interpolation::Histogram; 
    double x1 = 0.;
    double y1 = 4.;
    double x2 = 5.;
    double y2 = 8.;
    
    double x_low = x1;
    double x_hi = x2;
    double i = y1*(x_hi - x_low);
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = x2;
    x_hi = x1;
    EXPECT_DOUBLE_EQ(-i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = 1.3;
    x_hi = 4.7;
    i = y1*(x_hi - x_low);
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));
  }

  TEST(Integration, LinearLinear) {
    Interpolation interp = Interpolation::LinLin; 
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 5.;
    double y2 = 5.;
    
    double x_low = x1;
    double x_hi = x2;
    double i = 12.5;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = x2;
    x_hi = x1;
    EXPECT_DOUBLE_EQ(-i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = 1.3;
    x_hi = 4.7;
    i = 10.2;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));
  }

  TEST(Integration, LinearLog) {
    Interpolation interp = Interpolation::LinLog; 
    double x1 = 5.;
    double y1 = 1.;
    double x2 = 15.;
    double y2 = 21.;
    
    double x_low = x1;
    double x_hi = x2;
    double i = 127.95215467463257;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = x2;
    x_hi = x1;
    EXPECT_DOUBLE_EQ(-i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = 6.4;
    x_hi = 12.8;
    i = 80.169216888001344;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));
  }

  TEST(Integration, LogLinear) {
    Interpolation interp = Interpolation::LogLin; 
    double x1 = -10.;
    double y1 = 1.;
    double x2 = 15.;
    double y2 = 11.;
    
    double x_low = x1;
    double x_hi = x2;
    double i = 104.25809785606157;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = x2;
    x_hi = x1;
    EXPECT_DOUBLE_EQ(-i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = -4.378;
    x_hi = 5.69;
    i = 29.078502985642945;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));
  }

  TEST(Integration, LogLog) {
    Interpolation interp = Interpolation::LogLog; 
    double x1 = 1;
    double y1 = 5.;
    double x2 = 5.;
    double y2 = 8.;
    
    double x_low = x1;
    double x_hi = x2;
    double i = 27.089161107019226;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = x2;
    x_hi = x1;
    EXPECT_DOUBLE_EQ(-i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));

    x_low = 2.348;
    x_hi = 3.892;
    i = 10.739759942209645;
    EXPECT_DOUBLE_EQ(i, integrate(x_low, x_hi, x1, y1, x2, y2, interp));
  }

}
