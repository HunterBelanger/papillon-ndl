#include <PapillonNDL/region_1d.hpp>
#include <gtest/gtest.h>

namespace {
  using namespace pndl;

  TEST(Region1D, Min_Max_X) {
    std::vector<double> x {1.,2.,3.,4.,5.,6.};
    std::vector<double> y {1.,2.,3.,4.,5.,6.};
    Interpolation interp = Interpolation::LinLin;

    Region1D R(x,y,interp);

    EXPECT_DOUBLE_EQ(R.min_x(), x.front());
    EXPECT_DOUBLE_EQ(R.max_x(), x.back());
  }


  TEST(Region1D, Interpolation) {
    std::vector<double> x {1.,2.,3.,4.,5.,6.};
    std::vector<double> y {1.,2.,3.,4.,5.,6.};
    Interpolation interp = Interpolation::LogLin;

    Region1D R(x,y,interp);

    EXPECT_EQ(R.interpolation(), interp);
  }

  TEST(Region1D, EvaluationContinuous) {
    std::vector<double> x_vals {1.,2.,3.,4.,5.,6.};
    std::vector<double> y_vals {1.,2.,3.,4.,5.,6.};
    Interpolation interp = Interpolation::LinLin;

    Region1D R(x_vals,y_vals,interp);
    
    EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()-0.001));
    
    EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()));
    
    double x = 1.23456;
    EXPECT_DOUBLE_EQ(x, R(x));

    x = 3.238197263;
    EXPECT_DOUBLE_EQ(x, R(x));

    x = 4.238197263;
    EXPECT_DOUBLE_EQ(x, R(x));

    EXPECT_DOUBLE_EQ(x_vals.back(), R(x_vals.back()));

    EXPECT_DOUBLE_EQ(x_vals.back(), R(x_vals.back()+0.001));
  }

  TEST(Region1D, EvaluationDiscontinuous) {
    std::vector<double> x_vals {1.,2.,3.,3.,4.,5.};
    std::vector<double> y_vals {1.,2.,3.,6.,8.,10.};
    Interpolation interp = Interpolation::LinLin;

    Region1D R(x_vals,y_vals,interp);
    
    EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()-0.001));
    
    EXPECT_DOUBLE_EQ(x_vals.front(), R(x_vals.front()));
    
    double x = 1.23456;
    EXPECT_DOUBLE_EQ(x, R(x));

    x = 2.999999999;
    EXPECT_DOUBLE_EQ(x, R(x));

    x = 3.000000001;
    EXPECT_DOUBLE_EQ(2*x, R(x));

    x = 4.238197263;
    EXPECT_DOUBLE_EQ(2*x, R(x));

    EXPECT_DOUBLE_EQ(2*x_vals.back(), R(x_vals.back()));

    EXPECT_DOUBLE_EQ(2*x_vals.back(), R(x_vals.back()+0.001));
  }

  TEST(Region1D, Integration) {
    std::vector<double> x_vals {0.,2.,4.,6.};
    std::vector<double> y_vals {1.,1.,3.,2.};
    Interpolation interp = Interpolation::LinLin; 

    Region1D R(x_vals,y_vals,interp);

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
    EXPECT_DOUBLE_EQ(i1+i2+i3, R.integrate(x_low, x_hi));

    x_low = 3.;
    x_hi = 5.;
    double i = 5.25;
    EXPECT_DOUBLE_EQ(i, R.integrate(x_low, x_hi));
  }

  TEST(Region1D, Size) {
    std::vector<double> x_vals {1.,2.,3.,4.,5.,6.};
    std::vector<double> y_vals {1.,2.,3.,4.,5.,6.};
    Interpolation interp = Interpolation::LinLin;

    Region1D R(x_vals,y_vals,interp);

    EXPECT_EQ(x_vals.size(), R.size());
  }

  TEST(Region1D, X_Y) {
    std::vector<double> x_vals {1.,2.,3.,4.,5.,6.};
    std::vector<double> y_vals {1.,2.,3.,4.,5.,6.};
    Interpolation interp = Interpolation::LinLin;

    Region1D R(x_vals,y_vals,interp);

    const std::vector<double>& xref = R.x();
    const std::vector<double>& yref = R.y();

    EXPECT_EQ(xref.size(), yref.size());
    EXPECT_EQ(xref.size(), R.size());
    EXPECT_EQ(x_vals.size(), y_vals.size());
    EXPECT_EQ(x_vals.size(), xref.size());

    for(size_t i = 0; i < xref.size(); i++) {
      EXPECT_DOUBLE_EQ(xref[i], x_vals[i]);
      EXPECT_DOUBLE_EQ(yref[i], y_vals[i]); 
    }
  }

}
