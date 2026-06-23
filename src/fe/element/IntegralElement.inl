#include "fe/element/IntegralElement.h"

template <int num>
int EFEM::TetGaussIntegrationPoints<num>::getGaussPointNum()
{
    return num;
}

template <int num>
double EFEM::TetGaussIntegrationPoints<num>::getGaussPointWeight(int index)
{
    if (num == 1) {
        return 1.0;
    } else if (num == 4) {
        return 0.25;
    }
}

/*
 *                     3
 *                   ,/|`\
 *                 ,/  |  `\
 *               ,7    '.   `9
 *             ,/       8     `\
 *          ,/          |       `\
 *         0--------6--'.--------2
 *          `\.         |      ,/
 *             `\.      |    ,5
 *                `4.   '. ,/
 *                   `\. |/
 *                      `1
 */
template <int num>
Eigen::Vector3d EFEM::TetGaussIntegrationPoints<num>::getGaussPointLocalPosition(int index)
{
    Eigen::Vector3d shape_func;
    if (num == 1) {
        shape_func = Eigen::Vector3d(0.25, 0.25, 0.25); // 0.25, 0.25, 0.25, 0.25
    } else if (num == 4) {
        if (index == 0) {
            shape_func = Eigen::Vector3d(0.58541020, 0.13819660, 0.13819660); // 0.58541020, 0.13819660, 0.13819660, 0.13819660
        } else if (index == 1) {
            shape_func = Eigen::Vector3d(0.13819660, 0.58541020, 0.13819660); // 0.13819660, 0.58541020, 0.13819660, 0.13819660
        } else if (index == 2) {
            shape_func = Eigen::Vector3d(0.13819660, 0.13819660, 0.58541020); // 0.13819660, 0.13819660, 0.58541020, 0.13819660
        } else if (index == 3) {
            shape_func = Eigen::Vector3d(0.13819660, 0.13819660, 0.13819660); // 0.13819660, 0.13819660, 0.13819660, 0.58541020
        }
    }
    return shape_func;
}

template <int num>
int EFEM::HexGaussIntegrationPoints<num>::getGaussPointNum()
{
    return num;
}

template <int num>
double EFEM::HexGaussIntegrationPoints<num>::getGaussPointWeight(int index)
{
    double weight;
    if (num == 1) {
        return 1.0;
    } else if (num == 8) {
        weight = 0.125;
    } else if (num == 27) {
        if (index == 0) {
            weight = 125.00;
        } else if (index == 1) {
            weight = 200.00;
        } else if (index == 2) {
            weight = 125.00;
        } else if (index == 3) {
            weight = 200.00;
        } else if (index == 4) {
            weight = 320.00;
        } else if (index == 5) {
            weight = 200.00;
        } else if (index == 6) {
            weight = 125.00;
        } else if (index == 7) {
            weight = 200.00;
        } else if (index == 8) {
            weight = 125.00;
        } else if (index == 9) {
            weight = 200.00;
        } else if (index == 10) {
            weight = 320.00;
        } else if (index == 11) {
            weight = 200.00;
        } else if (index == 12) {
            weight = 320.00;
        } else if (index == 13) {
            weight = 512.00;
        } else if (index == 14) {
            weight = 320.00;
        } else if (index == 15) {
            weight = 200.00;
        } else if (index == 16) {
            weight = 320.00;
        } else if (index == 17) {
            weight = 200.00;
        } else if (index == 18) {
            weight = 125.00;
        } else if (index == 19) {
            weight = 200.00;
        } else if (index == 20) {
            weight = 125.00;
        } else if (index == 21) {
            weight = 200.00;
        } else if (index == 22) {
            weight = 320.00;
        } else if (index == 23) {
            weight = 200.00;
        } else if (index == 24) {
            weight = 125.00;
        } else if (index == 25) {
            weight = 200.00;
        } else if (index == 26) {
            weight = 125.00;
        }
        weight /= 5832.0;
    }
    return weight;
}

/*
 *      3----10----2
 *      |\         |\
 *      | 19       | 18
 *     11  \       9  \
 *      |   7----14+---6
 *      |   |      |   |
 *      0---+-8----1   |
 *       \  15      \  13
 *       16 |        17|
 *         \|         \|
 *          4----12----5
 */
template <int num>
Eigen::Vector3d EFEM::HexGaussIntegrationPoints<num>::getGaussPointLocalPosition(int index)
{
    Eigen::Vector3d shape_func;
    if (num == 1) {
        shape_func = Eigen::Vector3d(0.0, 0.0, 0.0);
    } else if (num == 8) {
        if (index == 0) {
            shape_func = Eigen::Vector3d(-1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0));
        } else if (index == 1) {
            shape_func = Eigen::Vector3d(1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0));
        } else if (index == 2) {
            shape_func = Eigen::Vector3d(1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0));
        } else if (index == 3) {
            shape_func = Eigen::Vector3d(-1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0));
        } else if (index == 4) {
            shape_func = Eigen::Vector3d(-1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0));
        } else if (index == 5) {
            shape_func = Eigen::Vector3d(1.00 / std::sqrt(3.0), -1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0));
        } else if (index == 6) {
            shape_func = Eigen::Vector3d(1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0));
        } else if (index == 7) {
            shape_func = Eigen::Vector3d(-1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0), 1.00 / std::sqrt(3.0));
        }
    } else if (num == 27) {
        if (index == 0) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 1) {
            shape_func = Eigen::Vector3d(0.0, -std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 2) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 3) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), 0.0, -std::sqrt(3.00 / 5.00));
        } else if (index == 4) {
            shape_func = Eigen::Vector3d(0.0, 0.0, -std::sqrt(3.00 / 5.00));
        } else if (index == 5) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), 0.0, -std::sqrt(3.00 / 5.00));
        } else if (index == 6) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 7) {
            shape_func = Eigen::Vector3d(0.0, std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 8) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00));
        } else if (index == 9) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 10) {
            shape_func = Eigen::Vector3d(0.0, -std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 11) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 12) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), 0.0, 0.0);
        } else if (index == 13) {
            shape_func = Eigen::Vector3d(0.0, 0.0, 0.0);
        } else if (index == 14) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), 0.0, 0.0);
        } else if (index == 15) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 16) {
            shape_func = Eigen::Vector3d(0.0, std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 17) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), 0.0);
        } else if (index == 18) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        } else if (index == 19) {
            shape_func = Eigen::Vector3d(0.0, -std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        } else if (index == 20) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), -std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        } else if (index == 21) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), 0.0, std::sqrt(3.00 / 5.00));
        } else if (index == 22) {
            shape_func = Eigen::Vector3d(0.0, 0.0, std::sqrt(3.00 / 5.00));
        } else if (index == 23) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), 0.0, std::sqrt(3.00 / 5.00));
        } else if (index == 24) {
            shape_func = Eigen::Vector3d(-std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        } else if (index == 25) {
            shape_func = Eigen::Vector3d(0.0, std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        } else if (index == 26) {
            shape_func = Eigen::Vector3d(std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00), std::sqrt(3.00 / 5.00));
        }
    }
    return shape_func;
}