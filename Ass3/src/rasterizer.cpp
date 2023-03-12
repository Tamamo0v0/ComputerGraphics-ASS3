#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    
    int len = sqrt(sample_rate);

    sample_buffer[y * width * len + x] = c;

  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);
    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    int len = sqrt(sample_rate);

    for(int i=0; i<len; ++i){
      for(int j=0; j<len; ++j){
        fill_pixel(sx * len + i, sy * len + j, color);
      }
    }

    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }
    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  bool inside_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    float px, float py) {

      float L0,L1,L2;
      float dX0, dY0, dX1, dY1, dX2, dY2;

      dX0 = x1 - x0;
      dY0 = y1 - y0;
      dX1 = x2 - x1;
      dY1 = y2 - y1;
      dX2 = x0 - x2;
      dY2 = y0 - y2;
      
      L0 = -(px - x0) * dY0 + (py - y0) * dX0;
      L1 = -(px - x1) * dY1 + (py - y1) * dX1;
      L2 = -(px - x2) * dY2 + (py - y2) * dX2;

      return ((L0 > 0) && (L1 > 0) && (L2 > 0)) || ((L0 < 0) && (L1 < 0) && (L2 < 0));

    }


  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // TODO: Task 2: Update to implement super-sampled rasterization
    int low_x, low_y, high_x, high_y;
    low_x = floor(min(x2, min(x0, x1)));
    low_y = floor(min(y2, min(y0, y1)));
    high_x = ceil(max(x2, max(x0, x1)));
    high_y = ceil(max(y2, max(y0, y1)));

    for(int i=low_x; i<=high_x; ++i){
      for(int j=low_y; j<=high_y; ++j){

        if (i < 0 || i >= width) continue;
        if (j < 0 || j >= height) continue;

        int len = sqrt(sample_rate);

        float interval = 1.0/(len + 1);

        for(int sub_x=1; sub_x<=len; ++sub_x){
          for(int sub_y=1; sub_y<=len; ++sub_y){

            float t_x,t_y;
            t_x = i + interval*sub_x;
            t_y = j + interval*sub_y;

            if(inside_triangle(x0,y0,x1,y1,x2,y2,t_x,t_y)){
              fill_pixel(i * len + sub_x - 1,j * len + sub_y - 1, color);
            }
          }
        }
    
      }
    }

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    int low_x, low_y, high_x, high_y;
    low_x = floor(min(x2, min(x0, x1)));
    low_y = floor(min(y2, min(y0, y1)));
    high_x = ceil(max(x2, max(x0, x1)));
    high_y = ceil(max(y2, max(y0, y1)));

    for(int i=low_x; i<=high_x; ++i){
      for(int j=low_y; j<=high_y; ++j){

        if (i < 0 || i >= width) continue;
        if (j < 0 || j >= height) continue;

        int len = sqrt(sample_rate);

        float interval = 1.0/(len + 1);

        for(int sub_x=1; sub_x<=len; ++sub_x){
          for(int sub_y=1; sub_y<=len; ++sub_y){
            float t_x,t_y;
            t_x = i + interval*sub_x;
            t_y = j + interval*sub_y;
            
            float u,v,w;
            
            Vector2D v0,v1,v2;

            v0 = {x1-x0,y1-y0};
            v1 = {x2-x0,y2-y0};
            v2 = {t_x-x0,t_y-y0};

            float d00 = dot(v0,v0);
            float d01 = dot(v0,v1);
            float d11 = dot(v1,v1);
            float d20 = dot(v2,v0);
            float d21 = dot(v2,v1);
            float denom = d00 * d11 - d01 * d01;
            
            v = (d11 * d20 - d01 * d21) / denom;
            w = (d00 * d21 - d01 * d20) / denom;
            u = 1.0 - v - w;

            if(v >= 0 && w >= 0 && (v + w) <= 1){
              Color color = u * c0 + v * c1 + w * c2;
              fill_pixel(i * len + sub_x - 1,j * len + sub_y - 1, color);
            }

          }
        }
    
      }
    }


  }

  Vector3D barycentric(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    float px, float py){

      float u,v,w;
            
      Vector2D vec0,vec1,vec2;

      vec0 = {x1-x0,y1-y0};
      vec1 = {x2-x0,y2-y0};
      vec2 = {px-x0,py-y0};

      float d00 = dot(vec0,vec0);
      float d01 = dot(vec0,vec1);
      float d11 = dot(vec1,vec1);
      float d20 = dot(vec2,vec0);
      float d21 = dot(vec2,vec1);
      float denom = d00 * d11 - d01 * d01;
      
      

      v = (d11 * d20 - d01 * d21) / denom;
      w = (d00 * d21 - d01 * d20) / denom;
      u = 1.0 - v - w;

      
      return {u,v,w};

    }



  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
  
    int low_x, low_y, high_x, high_y;
    low_x = floor(min(x2, min(x0, x1)));
    low_y = floor(min(y2, min(y0, y1)));
    high_x = ceil(max(x2, max(x0, x1)));
    high_y = ceil(max(y2, max(y0, y1)));
    Color col;

    for(int i=low_x; i<=high_x; ++i){
      for(int j=low_y; j<=high_y; ++j){

        if (i < 0 || i >= width) continue;
        if (j < 0 || j >= height) continue;

        int len = sqrt(sample_rate);

        float interval = 1.0/(len + 1);

        for(int sub_x=1; sub_x<=len; ++sub_x){
          for(int sub_y=1; sub_y<=len; ++sub_y){

            float t_x,t_y;
            t_x = i + interval*sub_x;
            t_y = j + interval*sub_y;
            
            float u,v,w;
            
            Vector3D uvw = barycentric(x0,y0,x1,y1,x2,y2,t_x,t_y);
            u = uvw[0];
            v = uvw[1];
            w = uvw[2];


            if(v >= 0 && w >= 0 && (v + w) <= 1){
              
              
              //cout << "u v are " << u0 << " " << v0 << endl;
              if(lsm == L_NEAREST || lsm == L_LINEAR){
                float u_dx,v_dx,w_dx;
                Vector3D uvw_dx = barycentric(x0,y0,x1,y1,x2,y2,t_x+1,t_y);
                u_dx = uvw_dx[0];
                v_dx = uvw_dx[1];
                w_dx = uvw_dx[2];


                float u_dy,v_dy,w_dy;
                Vector3D uvw_dy = barycentric(x0,y0,x1,y1,x2,y2,t_x,t_y+1);
                u_dy = uvw_dy[0];
                v_dy = uvw_dy[1];
                w_dy = uvw_dy[2];

                SampleParams param;
                param.psm = psm;
                param.lsm = lsm;

                param.p_uv = {u*u0 + v*u1 + w*u2, u*v0 + v*v1 + w*v2};
                param.p_dx_uv = {u_dx*u0 + v_dx*u1 + w_dx*u2, u_dx*v0 + v_dx*v1 + w_dx*v2};
                param.p_dy_uv = {u_dy*u0 + v_dy*u1 + w_dy*u2, u_dy*v0 + v_dy*v1 + w_dy*v2};
                col = tex.sample(param);
                fill_pixel(i * len + sub_x - 1,j * len + sub_y - 1, col);

              }else{
                float t_u,t_v;
                t_u = u * u0 + v * u1 + w * u2;
                t_v = u * v0 + v * v1 + w * v2;

                if(psm == P_NEAREST){
                  col = tex.sample_nearest({t_u,t_v});
                }else{
                  col = tex.sample_bilinear({t_u,t_v});
                }

                fill_pixel(i * len + sub_x - 1,j * len + sub_y - 1, col);
              }

            }
          }
        }
    
      }
    }

  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {

    unsigned int temp = floor(sqrt(rate));

    temp = temp * temp;

    if (temp != rate){
      cout << "The setting of sample rate isn't a square number, round it to the nearest one belows it."<<endl;
    }

    this->sample_rate = temp;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    int len = sqrt(sample_rate);

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = Color::Black;
      
        for(int i=0; i<len; ++i){
          for(int j=0; j<len; ++j){
            int temp_y = len * y + j;
            int temp_x = len * x + i;

            col += sample_buffer[width * len * temp_y + temp_x] * (1.0/sample_rate);
          }
        }

        //cout << "col = " << col.r << " " << col.g << " " << col.b << endl;

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }

      }
    }


  }

  Rasterizer::~Rasterizer() { }


}// CGL
