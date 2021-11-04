#pragma once

#include <glm/glm.hpp>

glm::mat4 inverse_orthogonal(const glm::mat4& m)
{
    // Based on compute_inverse<4, 4, T, Q, Aligned> from extern/glm/glm/detail/func_matrix.inl
    // a lot of values in orthogonal matrices (see extern/glm/glm/ext/matrix_clip_space.inl) are zero,
    // simplify compute_inverse<4, 4, T, Q, Aligned> using all known zeros of GLM_FUNC_QUALIFIER mat<4, 4, T, defaultp> ortho(T left, T right, T bottom, T top)

    auto Coef00 = m[2][2];
    auto Coef06 = m[1][1];

    auto Coef08 = -m[3][1] * m[2][2];
    auto Coef11 = m[1][1] * m[2][2];

    auto Coef16 = -m[3][0] * m[2][2];
    auto Coef22 = -m[3][0] * m[1][1];

    mat<4, 4, T, Q> Inverse(
        m[1][1]*Coef00, 0, 0, 0,
        0, m[0][0]*Coef00, 0, 0,
        0, 0, m[0][0]*Coef06, 0,
        m[1][1]*Coef16, m[0][0]*Coef08, 0, m[0][0]*Coef11
    );

    vec<4, T, Q> Row0(
        Inverse[0][0], 
        Inverse[1][0], 
        Inverse[2][0], 
        Inverse[3][0]
    );

    auto Determinant = (m[0][0] * Inverse[0][0]);
    T OneOverDeterminant = static_cast<T>(1) / Determinant;

    return Inverse * OneOverDeterminant;
}

