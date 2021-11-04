#pragma once

#include <glm/glm.hpp>

namespace inverse_ortho {
    // A lot of values in orthogonal matrices (see glm/glm/ext/matrix_clip_space.inl) are known.
    // Simplify compute_inverse<4, 4, T, Q, Aligned> (see glm/glm/detail/func_matrix.inl) using all known values of the respective orthogonal matrix.

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoLRBT(::glm::mat<4, 4, T, Q> const& m)
    {
        // glm/glm/ext/matrix_clip_space.inl
        // GLM_FUNC_QUALIFIER ::glm::mat<4, 4, T, defaultp> ortho(T left, T right, T bottom, T top)
        // ::glm::mat<4, 4, T, defaultp> Result(static_cast<T>(1));
        // Result[0][0] = static_cast<T>(2) / (right - left);
        // Result[1][1] = static_cast<T>(2) / (top - bottom);
        // Result[2][2] = - static_cast<T>(1);
        // Result[3][0] = - (right + left) / (right - left);
        // Result[3][1] = - (top + bottom) / (top - bottom);
        ::glm::mat<4, 4, T, Q> Inverse(
            -m[1][1], 0, 0, 0,
            0, -m[0][0], 0, 0,
            0, 0, m[0][0]*m[1][1], 0,
            m[1][1]*m[3][0], m[0][0]*m[3][1], 0, -m[0][0]*m[1][1]
        );
        
        T Determinant = (m[0][0] * Inverse[0][0]);
        T OneOverDeterminant = static_cast<T>(1) / Determinant;

        return Inverse * OneOverDeterminant;
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoLRBTNF(::glm::mat<4, 4, T, Q> const& m)
    {
        // glm/glm/ext/matrix_clip_space.inl
        // GLM_FUNC_QUALIFIER ::glm::mat<4, 4, T, defaultp> orthoLH_ZO(T left, T right, T bottom, T top, T zNear, T zFar)
        // ::glm::mat<4, 4, T, defaultp> Result(1);
        // Result[0][0] = static_cast<T>(2) / (right - left);
        // Result[1][1] = static_cast<T>(2) / (top - bottom);
        // Result[2][2] = static_cast<T>(1) / (zFar - zNear);
        // Result[3][0] = - (right + left) / (right - left);
        // Result[3][1] = - (top + bottom) / (top - bottom);
        // Result[3][2] = - zNear / (zFar - zNear);
        T Coef08 = -m[3][1] * m[2][2];
        T Coef10 = m[1][1] * m[3][2];
        T Coef11 = m[1][1] * m[2][2];
        T Coef16 = -m[3][0] * m[2][2];

        ::glm::mat<4, 4, T, Q> Inverse(
            m[1][1]*m[2][2], 0, 0, 0,
            0, m[0][0]*m[2][2], 0, 0,
            0, 0, m[0][0]*m[1][1],0,
            m[1][1]*Coef16, m[0][0]*Coef08, -m[0][0]*Coef10, m[0][0]*Coef11
        );

        T Determinant = (m[0][0] * Inverse[0][0]);
        T OneOverDeterminant = static_cast<T>(1) / Determinant;

        return Inverse * OneOverDeterminant;
    }

    // remaining functions have the same implementation as inverse_orthoLRBTNF:

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoLH_ZO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoLH_NO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoRH_ZO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoRH_NO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoZO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoNO(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoLH(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_orthoRH(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

    template<typename T, ::glm::qualifier Q>
    GLM_FUNC_DECL ::glm::mat<4, 4, T, Q> inverse_ortho(::glm::mat<4, 4, T, Q> const& m)
    {
        return inverse_orthoLRBTNF(m);
    }

} // namespace inverse_ortho

