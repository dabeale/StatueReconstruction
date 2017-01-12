#include "libppm.h"

namespace PPM
{
    void Write( const char* filename, const uint32_t height, const uint32_t width, const double* image )
    {
        std::ofstream ppm_out(filename, std::ios::binary);

        ppm_out<<"P6";
        ppm_out<<' ';
        ppm_out<<width;
        ppm_out<<' ';
        ppm_out<<height;
        ppm_out<<' ';
        ppm_out<<"255";
        ppm_out<<std::endl;

        for(unsigned y=0;y<height;y++)
        for(unsigned x=0;x<width;x++)
        {
          const double d = image[width*y + x]* 255.0;

          const char r = (d>254) ? 254 : ((d<0)?0:static_cast<uint32_t>(d));
          const char g = r;
          const char b = r;

          ppm_out<<r<<g<<b;
        }

        ppm_out.flush();
        ppm_out.close();
    }

    std::vector<double> Read(const char* filename, uint32_t &width, uint32_t &height)
    {
        std::ifstream ppm_in(filename,std::ios::binary);

        std::string magic_number("  ");

        ppm_in.get(magic_number[0]);
        ppm_in.get(magic_number[1]);

        if (magic_number != std::string("P6"))
        {
            std::cerr<<"error: unrecognized file format\n"<<filename<<" is not a PPM file.\n"<<std::endl;
            exit(2);
        }

        unsigned bpp;

        ppm_in>>width>>height>>bpp;
        std::vector<double> out(width*height);
        double* image = out.data();

        if (bpp != 255)
        {
            std::cerr<<"error: unsupported maximum value ("<<bpp<<")\n"<<"It must be 255."<<std::endl;
            exit(3);
        }

        char ch;
        ppm_in.get(ch); // Trailing white space.

        char r,g,b;

        for(unsigned y=0; y<height; y++)
        for(unsigned x=0; x<width; x++)
        {
            ppm_in.get(r);
            ppm_in.get(g);
            ppm_in.get(b);

            const unsigned char R = static_cast<unsigned char>(r);
            const unsigned char G = static_cast<unsigned char>(g);
            const unsigned char B = static_cast<unsigned char>(b);

            image[y*width +x] = (20.0 * R + 40.0 * G + 1.0 * B) / (61.0 * 255.0 );
        }

        return out;
    }
}
