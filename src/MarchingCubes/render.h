#ifndef RENDER_H
#define RENDER_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#include <vector>
#include <functional>
#include "marchingcubes.h"

class Render
{
public:
    struct MCParams
    {
        double isoval;
        uint32_t N;
        double max[3];
        double min[3];
    };

    enum MergeType
    {
        Add,
        Union,
        Bernhardt,
        Barthe,
        Rockwood,
        Functional,
        nvals
    };

    /**
     * @brief ImplicitFunction
     * An implicit function of a three vector and time
     */
    typedef std::function< double(const double*, const double) > ImplicitFunction;

    Render();
    Render(const MCParams mcp);
    ~Render();

    void AddFunction( const ImplicitFunction ifunc );
    void IncrementMergeType();
    void IncrementTime(const double dt);
    void SwitchSpin();

    void March();

    float GetPitch() const;
    float GetYaw() const;

    uint32_t GetNumberOfFunctions() const;

private:
    std::vector<ImplicitFunction> m_implicitFuncs; ///< A container of implicit functions
    ImplicitFunction m_renderFunction;

    MergeType m_mt; ///< The type of implicit merging to do
    MC::MarchingCubes m_mc; ///< A marching cubes object
    uint32_t m_currentType;
    bool m_spin;
    double m_t; ///< The current time

    float m_pitch, m_yaw;


    void SwitchFunc(MergeType mt);
};

namespace GLR
{
    void Resize( GLsizei iWidth, GLsizei iHeight ); ///< Resize the window
    void Keyboard(unsigned char cKey, int iX, int iY);
    void Special(int iKey, int iX, int iY);
    void Idle();
    void DrawScene();

    void Create(const Render::MCParams mcp);
    void AddFunction(const Render::ImplicitFunction ifunc );
    void Start();

    static Render::MCParams m_mcp;
    static Render* renderptr = NULL;
    static Render& GetRender()
    {
         if(renderptr == NULL)
         {
            renderptr = new Render(m_mcp);
         }
         return *renderptr;
    }
}

#endif // RENDER_H
