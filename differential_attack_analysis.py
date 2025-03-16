import numpy as np

def npcr(image1, image2):
  if image1.shape != image2.shape:
    raise ValueError("Images must be the same size.")
  num_changed_pixels = np.sum(image1 != image2)
  npcr = num_changed_pixels / (image1.shape[0] * image1.shape[1]) * 100

  return npcr


def uaci(img1,img2):
    h,w = img1.shape
    val=0.0
    for y in range(h):
        for x in range(w):
            val += abs(int(img1[x,y] - int(img2[x,y])))
    val=val/(w*h*255)*100
    return val
   




