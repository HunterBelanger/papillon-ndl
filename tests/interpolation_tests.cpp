#include <PapillonNDL/interpolation.hpp>
#include <gtest/gtest.h>

namespace {
  using namespace pndl;

  TEST(Interpolation, Histogram) {
    Interpolation interp = Interpolation::Histogram; 

    double x1 = 8.;
    double y1 = 4.;

    double x2 = 10.;
    double y2 = 6.;

    EXPECT_DOUBLE_EQ(y1, interpolate(x1, x1, y1, x2, y2, interp));
    
    double x = 9.;
    double y = y1;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));

    x = x2 - 0.000000001;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));

    EXPECT_DOUBLE_EQ(y2, interpolate(x2, x1, y1, x2, y2, interp));
  }

  TEST(Interpolation, LinearLinear) {
    Interpolation interp = Interpolation::LinLin;
    // Here, y is linear in x
    // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1

    double x1 = 8.;
    double y1 = 0.;

    double x2 = 10.;
    double y2 = 6.;
    EXPECT_DOUBLE_EQ(y1, interpolate(x1, x1, y1, x2, y2, interp));
    
    double x = 9.;
    double y = 3.;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));
    
    EXPECT_DOUBLE_EQ(y2, interpolate(x2, x1, y1, x2, y2, interp));
  }

  TEST(Interpolation, LinearLog) {
    // Here, y is linear in log(x)
    // y = (log(x/x1)/log(x2/x1))*(y2 - y1) + y1
    Interpolation interp = Interpolation::LinLog;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interpolate(x1, x1, y1, x2, y2, interp));
    
    double x = 2.;
    double y = 2.;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));

    x = 3.;
    y = 2.584962500721156;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));
    
    EXPECT_DOUBLE_EQ(y2, interpolate(x2, x1, y1, x2, y2, interp));
  }

  TEST(Interpolation, LogLinear) {
    // Here, y is linear in log(x)
    // log(y) = ((x - x1)/(x2 - x1))*log(y2/y1) + log(y1)
    Interpolation interp = Interpolation::LogLin;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interpolate(x1, x1, y1, x2, y2, interp));

    double x = 2.;
    double y = 1.4422495703074083;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));

    x = 3.;
    y = 2.080083823051904;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));
    
    
    EXPECT_DOUBLE_EQ(y2, interpolate(x2, x1, y1, x2, y2, interp));
  }

  TEST(Interpolation, LogLog) {
    // Here, log(y) is linear in log(x)
    // log(y) = (log(x/x1)/log(x2/x1))*log(y2/y1) + log(y1)
    Interpolation interp = Interpolation::LogLog;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interpolate(x1, x1, y1, x2, y2, interp));

    double x = 2.;
    double y =1.7320508075688774 ;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));

    x = 3.;
    y = 2.388414221757005;
    EXPECT_DOUBLE_EQ(y, interpolate(x, x1, y1, x2, y2, interp));
    
    
    EXPECT_DOUBLE_EQ(y2, interpolate(x2, x1, y1, x2, y2, interp));
  }

}
