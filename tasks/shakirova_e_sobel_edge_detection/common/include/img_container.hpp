#pragma once

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <stb/stb_image.h>

namespace shakirova_e_sobel_edge_detection {

struct ImgContainer {
  int width{0};
  int height{0};
  std::vector<int> pixels;

  ImgContainer() = default;
  ImgContainer(int w, int h) : width(w), height(h), pixels(static_cast<size_t>(w) * static_cast<size_t>(h), 0) {}
  ImgContainer(int w, int h, std::vector<int> px)
      : width(w), height(h), pixels(std::move(px)) {}

  [[nodiscard]] bool IsValid() const {
    return width > 2 && height > 2 &&
           static_cast<int>(pixels.size()) == width * height;
  }

  [[nodiscard]] int& At(int x, int y) { return pixels[(y * width) + x]; }
  [[nodiscard]] const int& At(int x, int y) const { return pixels[(y * width) + x]; }

  friend bool operator==(const ImgContainer& lhs, const ImgContainer& rhs) {
    return lhs.width == rhs.width && lhs.height == rhs.height &&
           lhs.pixels == rhs.pixels;
  }

  static ImgContainer FromFile(const std::string& path) {
    if (!std::filesystem::exists(path)) {
      throw std::runtime_error("File not found: " + path);
    }
    const std::string ext = GetExt(path);
    if (ext == ".png" || ext == ".jpg" || ext == ".jpeg" || ext == ".bmp") {
      return LoadRaster(path);
    }
    if (ext == ".txt") {
      return LoadTXT(path);
    }
    throw std::runtime_error("Unsupported format: " + ext +
                              " (supported: .png .jpg .bmp .txt)");
  }

 private:
  static int ToGray(uint8_t r, uint8_t g, uint8_t b) {
    return static_cast<int>((0.299 * r) + (0.587 * g) + (0.114 * b));
  }

  static std::string GetExt(const std::string& path) {
    const auto pos = path.rfind('.');
    if (pos == std::string::npos) {
      return "";
    }
    std::string ext = path.substr(pos);
    std::ranges::transform(ext, ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext;
  }

  static ImgContainer LoadRaster(const std::string& path) {
    int w = 0;
    int h = 0;
    int ch = 0;
    uint8_t* data = stbi_load(path.c_str(), &w, &h, &ch, STBI_rgb);
    if (data == nullptr) {
      throw std::runtime_error("stb_image failed to load: " + path);
    }
    ImgContainer img(w, h);
    for (int i = 0; i < w * h; ++i) {
      img.pixels[i] = ToGray(data[static_cast<ptrdiff_t>(i) * 3],
                              data[(static_cast<ptrdiff_t>(i) * 3) + 1],
                              data[(static_cast<ptrdiff_t>(i) * 3) + 2]);
    }
    stbi_image_free(data);
    return img;
  }

  static ImgContainer LoadTXT(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
      throw std::runtime_error("Cannot open file: " + path);
    }
    int w = 0;
    int h = 0;
    f >> w >> h;
    if (w <= 0 || h <= 0) {
      throw std::runtime_error("Invalid dimensions in: " + path);
    }
    ImgContainer img(w, h);
    for (int i = 0; i < w * h; ++i) {
      f >> img.pixels[i];
    }
    return img;
  }
};

}  // namespace shakirova_e_sobel_edge_detection