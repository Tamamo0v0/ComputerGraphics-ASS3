#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    float level = get_level(sp);
    int level_rounded = (int)round(level);
    //cout << "sample" << endl;
    if(sp.lsm == L_NEAREST){
      //cout << "L_near" << endl;
      if (sp.psm == P_NEAREST)
      {
        //cout << "P_near" << endl;
        return sample_nearest(sp.p_uv, level_rounded);
      }
      else if (sp.psm == P_LINEAR)
      {
        //cout << "P_linear" << endl;
        return sample_bilinear(sp.p_uv, level_rounded);
      }
    }else if(sp.lsm == L_LINEAR){
      //cout << "L_liear" << endl;
      if (sp.psm == P_NEAREST)
      {
        //cout << "P_near" << endl;
        float low_weight = level - floor(level);
        float high_weight = ceil(level) - level;
        Color low_c, high_c, col;
        low_c = sample_nearest(sp.p_uv, (int)floor(level));
        high_c = sample_nearest(sp.p_uv, (int)ceil(level));
        col = low_weight * low_c + high_weight * high_c;

        return col;
      }
      else if (sp.psm == P_LINEAR)
      {
        //cout << "P_linear" << endl;
        float low_weight = level - floor(level);
        float high_weight = ceil(level) - level;
        Color low_c, high_c, col;
        low_c = sample_bilinear(sp.p_uv, (int)floor(level));
        high_c = sample_bilinear(sp.p_uv, (int)ceil(level));
        col = low_weight * low_c + high_weight * high_c;
        
        return col;
      }
    }
    


// return magenta for invalid level
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    auto& mip = mipmap[0];

    Vector2D p_uv = sp.p_uv;
    Vector2D p_dx_uv = sp.p_dx_uv;
    Vector2D p_dy_uv = sp.p_dy_uv;

    float dx,dy;
    //cout << "dx_uv" <<p_dx_uv[0] <<endl;
    dx = pow(pow((p_dx_uv[0] - p_uv[0]) * mip.width, 2) + pow((p_dx_uv[1] - p_uv[1]) * mip.height, 2), 0.5);
    dy = pow(pow((p_dy_uv[0] - p_uv[0]) * mip.width, 2) + pow((p_dy_uv[1] - p_uv[1]) * mip.height, 2), 0.5);

    float level = log(max(dx,dy))/log(2);

    //cout << "dx" << dx << "dy" << dy << "level" << level << endl;
    return level;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    if (level == -1) level = 0;
    //cout << "1" << endl;
    auto& mip = mipmap[level];
    //cout << "2" << endl;
    Color col;
    //cout << level << endl;
    int u = (int)round(uv[0] * mip.width);
    //cout << "1" << endl;
    int v = (int)round(uv[1] * mip.height);

    //cout << "2.5" << endl;
    // check bounds
    //if (u < 0 || u >= mip.width) return Color(1, 0, 1);
    //if (v < 0 || v >= mip.height) return Color(1, 0, 1);

    //cout << "3" << endl;
    col = mip.get_texel(u,v);
    //cout << "4" << endl;
    return col;
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    if (level < 0) level = 0;
    auto& mip = mipmap[level];

    Color col;
    int u = (int)round(uv[0] * mip.width);
    int v = (int)round(uv[1] * mip.height);
    // check bounds
    //if (u < 0 || u >= mip.width) return Color(1, 0, 1);
    //if (v < 0 || v >= mip.height) return Color(1, 0, 1);

    Color u00,u01,u10,u11;
    Color u0,u1;

    u00 = mip.get_texel(u,v);
    u10 = mip.get_texel(u+1,v);
    u01 = mip.get_texel(u,v+1);
    u11 = mip.get_texel(u+1,v+1);

    u0 = abs(uv[0] * mip.width - u) * u00 + abs(1 - abs(uv[0] * mip.width - u))* u10;
    u1 = abs(uv[0] * mip.width - u) * u01 + abs(1 - abs(uv[0] * mip.width - u)) * u11;

    col = abs(uv[1] * mip.height - v) * u0 + abs(1 - abs(uv[1] * mip.height - v)) * u1;

    return col;
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
