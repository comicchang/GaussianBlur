#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CHANNELS 3

typedef unsigned char uchar;

uchar *bilinearInterpolation(const uchar *image, int width, int height,
                             int newWidth, int newHeight);
void processImage(const char *inputPath, const char *outputPath, float radius);
void gaussianBlur(uchar *image, int width, int height, double sigma);
void convolutionX(const uchar *inputImage, uchar *outputImage, int width,
                  int height, const float *kernel, int kernelSize);
void convolutionY(const uchar *inputImage, uchar *outputImage, int width,
                  int height, const float *kernel, int kernelSize);

int calcOffsetWithClamp(int width, int height, int x, int y) {
  // Clamp x and y values to ensure they are within image bounds
  x = fmax(0, fmin(x, width - 1));
  y = fmax(0, fmin(y, height - 1));

  // Calculate the offset of the pixel
  int offset = (y * width + x) * CHANNELS;

  // Return the offset
  return offset;
}

int calcOffset(int width, int height, int x, int y) {
  // Calculate the offset of the pixel
  return (y * width + x) * CHANNELS;
}

uchar *bilinearInterpolation(const uchar *image, int width, int height,
                             int newWidth, int newHeight) {
  int channels = CHANNELS;
  float xRatio = (float)width / newWidth;
  float yRatio = (float)height / newHeight;

  uchar *outputImage =
      (uchar *)malloc(newWidth * newHeight * channels * sizeof(uchar));

  for (int newY = 0; newY < newHeight; newY++) {
    for (int newX = 0; newX < newWidth; newX++) {
      float originalX = newX * xRatio;
      float originalY = newY * yRatio;
      int x1 = (int)originalX;
      int y1 = (int)originalY;
      float alpha = originalX - x1;
      float beta = originalY - y1;

      int leftTopOffset = calcOffsetWithClamp(width, height, x1, y1);
      int outputOffset = calcOffset(newWidth, newHeight, newX, newY);

      for (int channel = 0; channel < channels; channel++) {
        int offset = leftTopOffset + channel;

        uchar leftTop = image[offset];
        uchar rightTop = image[offset + channels];
        uchar leftBottom = image[offset + width * channels];
        uchar rightBottom = image[offset + (width + 1) * channels];

        float interpolatedValue =
            (1 - alpha) * (1 - beta) * leftTop + alpha * (1 - beta) * rightTop +
            (1 - alpha) * beta * leftBottom + alpha * beta * rightBottom;

        outputImage[outputOffset + channel] = round(interpolatedValue);
      }
    }
  }

  return outputImage;
}

void fill1DGaussianKernel(float *kernel, int size, float sigma) {
  // 计算高斯核的中心位置
  int center = size / 2;

  // 计算高斯分布的标准差的平方
  float sigma2 = sigma * sigma;

  // 计算高斯核的总和，用于归一化
  float sum = 0.0f;

  // 填充高斯核
  for (int i = 0; i < size; i++) {
    // 计算当前位置与中心位置的差值
    int offset = i - center;

    // 计算高斯核的值
    float value = exp(-(offset * offset) / (2 * sigma2));

    // 将值存储到高斯核数组中
    kernel[i] = value;

    // 累加高斯核的值
    sum += value;
  }

  // 归一化高斯核
  for (int i = 0; i < size; i++) {
    kernel[i] /= sum;
  }

  // printf("Gaussian kernel: ");
  // for (int i = 0; i < kernelSize; i++) {
  //   printf("%f ", kernel[i]);
  // }
  // printf("\n");
}

void convolutionY(const uchar *inputImage, uchar *outputImage, int width,
                  int height, const float *kernel, int kernelSize) {
  int channels = CHANNELS;
  int kernelRadius = kernelSize / 2;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      float sum[channels];
      memset(sum, 0, sizeof(float) * channels);

      for (int i = 0; i < kernelSize; i++) {
        int offsetY = i - kernelRadius;
        int inputOffset = calcOffsetWithClamp(width, height, x, y + offsetY);

        for (int channel = 0; channel < channels; channel++) {
          sum[channel] += inputImage[inputOffset + channel] * kernel[i];
        }
      }

      uchar *output = outputImage + calcOffset(width, height, x, y);
      for (int channel = 0; channel < channels; channel++) {
        output[channel] = round(sum[channel]);
      }
    }
  }
}

void convolutionX(const uchar *inputImage, uchar *outputImage, int width,
                  int height, const float *kernel, int kernelSize) {
  int channels = CHANNELS;
  int kernelRadius = kernelSize / 2;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      float sum[channels];
      for (int channel = 0; channel < channels; channel++) {
        sum[channel] = 0;
      }

      for (int i = 0; i < kernelSize; i++) {
        int offsetX = i - kernelRadius;
        int inputOffset = calcOffsetWithClamp(width, height, x + offsetX, y);

        for (int channel = 0; channel < channels; channel++) {
          sum[channel] += inputImage[inputOffset + channel] * kernel[i];
        }
      }

      uchar *output = outputImage + calcOffset(width, height, x, y);
      for (int channel = 0; channel < channels; channel++) {
        output[channel] = round(sum[channel]);
      }
    }
  }
}

// void meanColor(uchar *image, int width, int height) {
//   int channels = CHANNELS;
//   float meanColor[CHANNELS];
//   for (int channel = 0; channel < CHANNELS; channel++) {
//     meanColor[channel] = 0;
//   }
//   for (int y = 0; y < height; y++) {
//     for (int x = 0; x < width; x++) {
//       for (int channel = 0; channel < channels; channel++) {
//         meanColor[channel] += image[calcOffset(width, height, x, y) + channel];
//       }
//     }
//   }
//   for (int channel = 0; channel < channels; channel++) {
//     printf("meanColor[%d] = %f\n", channel,
//            meanColor[channel] / width / height);
//   }
// }

void gaussianBlur(uchar *image, int width, int height, double sigma) {
  // Calculate kernel size based on sigma
  int kernelSize = 31;

  printf("Kernel size: %d\n", kernelSize);
  // Allocate memory for the kernel
  float *kernel = (float *)malloc(kernelSize * sizeof(float));

  // Fill the 1D Gaussian kernel
  fill1DGaussianKernel(kernel, kernelSize, sigma);

  uchar *temp = (uchar *)malloc(width * height * CHANNELS * sizeof(uchar));

  // meanColor(image, width, height);

  printf("ConvolutionX...size[%d, %d], kernel size %d\n", width, height,
         kernelSize);
  // Apply Gaussian blur horizontally
  convolutionX(image, temp, width, height, kernel, kernelSize);

  // meanColor(temp, width, height);

  printf("ConvolutionY...size[%d, %d], kernel size %d\n", width, height,
         kernelSize);
  // Apply Gaussian blur vertically
  convolutionY(temp, image, width, height, kernel, kernelSize);

  // meanColor(image, width, height);

  printf("Convolution...Done\n");

  // Free memory
  free(kernel);
  free(temp);
}

void processImage(const char *inputPath, const char *outputPath, float sigma) {
  // Load image using stb_image
  int origWidth, origHeight, channels;
  uchar *image =
      stbi_load(inputPath, &origWidth, &origHeight, &channels, CHANNELS);
  if (image == NULL) {
    printf("Failed to load image: %s\n", inputPath);
    return;
  }
  printf("Loaded image: %s, image size: %dx%d, channels: %d\n", inputPath,
         origWidth, origHeight, channels);

  // Apply Gaussian blur iteratively
  int width = origWidth;
  int height = origHeight;
  while (sigma > 4.0 && width >= 2 && height >= 2) {
    // Resize image using bilinear interpolation
    int newWidth = width / 2;
    int newHeight = height / 2;
    printf("Resizing image from %dx%d to %dx%d\n", width, height, newWidth,
           newHeight);
    uchar *resizedImage =
        bilinearInterpolation(image, width, height, newWidth, newHeight);
    free(image);
    image = resizedImage;
    width = newWidth;
    height = newHeight;
    sigma = sigma / 2.0;

    // Save the resized image
    // filename = "output + sigma + .png"
    char filename[1024];
    snprintf(filename, 1024, "%s_%g.png", outputPath, sigma);
    stbi_write_png(filename, width, height, CHANNELS, image, width * CHANNELS);
  }

  // Apply Gaussian blur
  gaussianBlur(image, width, height, sigma);

  char filename[1024];
  snprintf(filename, 1024, "%s_%g_blurred.png", outputPath, sigma);
  stbi_write_png(filename, width, height, CHANNELS, image, width * CHANNELS);

  printf("Resizing image from %dx%d to %dx%d\n", width, height, origWidth,
         origHeight);
  // Upscale the blurred image
  uchar *outputImage =
      bilinearInterpolation(image, width, height, origWidth, origHeight);
  free(image);

  printf("Saving image to %s\n", outputPath);
  // Save the final result
  stbi_write_png(outputPath, origWidth, origHeight, CHANNELS, outputImage,
                 origWidth * CHANNELS);

  // Free memory
  free(outputImage);
}

int main() {
  const char *inputPath = "../input.png";
  const char *outputPath = "output.png";
  float sigma = 30.0f;

  processImage(inputPath, outputPath, sigma);

  return 0;
}