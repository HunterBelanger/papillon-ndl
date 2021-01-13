#include <PapillonNDL/multi_region_1d.hpp>
#include <gtest/gtest.h>

#include <iostream>

namespace {
  using namespace pndl;

  TEST(MultiRegion1D, ConstructorRegions) {
    std::vector<std::shared_ptr<Region1D>> regions 
      {build_Region1D({1.,4.},{1., 2.},Interpolation::LinLin),
       build_Region1D({5.,6.},{2., 23.},Interpolation::LinLin)}; 
    EXPECT_ANY_THROW(MultiRegion1D m_fail(regions));

    std::vector<std::shared_ptr<Region1D>> regions2 
      {build_Region1D({1.,4.},{1., 2.},Interpolation::LinLin),
       build_Region1D({-1.,6.},{2., 23.}, Interpolation::LinLin)}; 
    EXPECT_ANY_THROW(MultiRegion1D m_fail2(regions2));

    std::vector<std::shared_ptr<Region1D>> regions3
      {build_Region1D({1.,4.},{1., 2.}, Interpolation::LinLin),
       build_Region1D({4.,6.},{2., 23.},Interpolation::LinLin)}; 
    EXPECT_NO_THROW(MultiRegion1D m(regions3));
  }

  TEST(MultiRegion1D, ConstructorACE) {
    std::vector<uint32_t> NBT1 {2,4,5};
    std::vector<Interpolation> INT1 {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x1 {1.,2.,6.,10.};
    std::vector<double> y1 {0.,1.,2.,6.,20.};
    EXPECT_ANY_THROW(MultiRegion1D R(NBT1, INT1, x1, y1));

    std::vector<uint32_t> NBT2 {2,5};
    std::vector<Interpolation> INT2 {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x2 {1.,2.,2.,6.,10.};
    std::vector<double> y2 {0.,1.,2.,6.,20.};
    EXPECT_ANY_THROW(MultiRegion1D R(NBT2, INT2, x2, y2));

    std::vector<uint32_t> NBT3 {2,4,5};
    std::vector<Interpolation> INT3 {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x3 {1.,2.,2.,1.5,10.};
    std::vector<double> y3 {0.,1.,2.,6.,20.};
    EXPECT_ANY_THROW(MultiRegion1D R(NBT3, INT3, x3, y3));

    std::vector<uint32_t> NBT {2,4,5};
    std::vector<Interpolation> INT {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x {1.,2.,2.,6.,10.};
    std::vector<double> y {0.,1.,2.,6.,20.};
    EXPECT_NO_THROW(MultiRegion1D R(NBT, INT, x, y));
  }

  TEST(MultiRegion1D, Size) {
    std::vector<std::shared_ptr<Region1D>> regions 
      {build_Region1D({1.,4.},{1., 2.},Interpolation::LinLin),
       build_Region1D({4.,6.},{2., 23.},Interpolation::LinLin)}; 
    MultiRegion1D R(regions);

    EXPECT_EQ(R.size(), regions.size());
  }

  TEST(MultiRegion1D, Region) {
    std::vector<std::shared_ptr<Region1D>> regions 
      {build_Region1D({1.,4.},{1., 2.},Interpolation::LinLin),
       build_Region1D({4.,6.},{2., 23.},Interpolation::LinLin)}; 
    MultiRegion1D R(regions);

    EXPECT_EQ(R.size(), regions.size());
  }

  TEST(MultiRegion1D, Min_Max_x) {
    std::vector<std::shared_ptr<Region1D>> regions 
      {build_Region1D({1.,4.},{1., 2.},Interpolation::LinLin),
       build_Region1D({4.,6.},{2., 23.},Interpolation::LinLin)}; 
    MultiRegion1D R(regions);

    EXPECT_DOUBLE_EQ(regions.front()->min_x(), R.min_x());
    EXPECT_DOUBLE_EQ(regions.back()->max_x(), R.max_x());
  }

  TEST(MultiRegion1D, Evaluation) {
    std::vector<std::shared_ptr<Region1D>> regions 
      {build_Region1D({1.,4.},{1., 4.},Interpolation::LinLin),
       build_Region1D({4.,8.},{8., 16.},Interpolation::LinLin)}; 
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
    y = 2*x;
    EXPECT_DOUBLE_EQ(y, R(x));

    x = 7.;
    y = 2*x;
    EXPECT_DOUBLE_EQ(y, R(x));

    x = 8.;
    y = 2*x;
    EXPECT_DOUBLE_EQ(y, R(x));

    x = 12.;
    EXPECT_DOUBLE_EQ(y, R(x));
  }

  TEST(MultiRegion1D, EvaluationDiscontinuousMuliregion) {
    std::vector<uint32_t> NBT {2,4,5};
    std::vector<Interpolation> INT {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x_vals {1.,2.,2.,6.,10.};
    std::vector<double> y_vals {0.,1.,2.,6.,20.}; 
    MultiRegion1D R(NBT,INT,x_vals,y_vals);
    
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
    std::vector<uint32_t> NBT {2,3,4};
    std::vector<Interpolation> INT {Interpolation::LinLin, Interpolation::LinLin, Interpolation::LinLin};
    std::vector<double> x_vals {0.,2.,4.,6.};
    std::vector<double> y_vals {1.,1.,3.,2.};

    MultiRegion1D R(NBT, INT, x_vals,y_vals);
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
    EXPECT_DOUBLE_EQ(i1+i2+i3, R.integrate(x_low, x_hi));

    x_low = 3.;
    x_hi = 5.;
    double i = 5.25;
    EXPECT_DOUBLE_EQ(i, R.integrate(x_low, x_hi));
  }

}
