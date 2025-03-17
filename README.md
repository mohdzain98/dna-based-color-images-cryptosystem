# Design of DNA based Color Images Cryptosystem

## Introduction

DNA based color images refer to the encoding of the color images into DNA sequences called as DNA computing.  
The basic idea behind DNA computing is that DNA can be used to store and process information in a way that is similar to how a computer uses binary digits (bits) to represent data.
Instead of representing images using traditional pixel values the image data is converted into a sequence of DNA nucleotides(Adenine-A, Thymine-T, Cytosine-C and Guanine-G).
For each color channel i.e. RGB of pixel might be represented by a certain DNA bases or combination of bases.
<strong>Overview of Conversion</strong>- The pixel value 0-255 can be converted into DNA sequence. Using a set of predefined rules we can map every two bit of a 8 bit binary sequence to one of the DNA nucleotides.

## The Proposed Cryptosystem

The cryptosystem makes use of six 2D-Multiple Collapse Chaotic Maps, which are produced by the Key Scheduling Algorithm from a 240-symmetric key. The first stage involves intermixing of channels of the image then utilising the proposed operator ⊗ to operate the rows of the intermixed images with the initial array. Using a modified version of Arnold's Cat Map [10], which can be applied to both square and rectangular images, the pixels of the composite image are jumbled. The suggested operator is then used to further mix the pixels in the directions of forward spiral row, forward spiral column, reverse spiral row, and reverse spiral column. These keyless procedures are necessary to eliminate the link between pixels and the image's visual perception. Lastly, a non-linear layer has been added by combining the suggested operator with the idea of DNA cryptography. One of the eight encoding principle identified by the rule chart is applied to the mixed image that was obtained in the preceding phases, converting it into a sequence of DNA nucleotides. At this stage DNA encoded image is obtained and these encodings are further substituted using proposed operator with the substitution map. Lastly, one of the eight DNA decoding guidelines identified by a different rule chart is applied to decode the substituted DNA sequences.

## Encryption Algorithm Flow Chart

<img src="images/efchart.jpg" alt="echart"/>

## Encryption Process

1. Inter Channel Mixing
2. Mix rows with Initial Vector
3. Arnold Cat Map algorithm
4. Spiral Mixing in four different directions
5. DNA Encoding
6. DNA Substitution
7. DNA Decoding
8. Final cipher image

## Decryption

Reverse of Encryption procedure.

## Experimental Results and Analysis

### Metrics

1. Key Space and Key Sensitivity.
2. Histogram Analysis
3. Variance Analysis
4. chi Square Analysis
5. Entropy Analysis
6. Correlation Analysis
7. NPCR and UACI
8. PSNR
9. Comparative Analysis with other operator and Referenced Present Cryptosystems.

## Set up and Getting Started

1. Clone the Repository:

   ```bash
   git clone https://github.com/mohdzain98/dna-based-color-images-cryptosystem.git
   cd dna-based-color-images-cryptosystem

   ```

2. Create Virtual Environment

   ```bash
   python -m veny venv_name

   ```

3. Install Dependencies

   - Windows
     ```bash
     venv/Scripts/activate
     ```
   - Linux/Mac
     ```bash
     source venv/bin/activate
     ```

   ```bash
   pip install -r requirements.txt

   ```

4. Running the cryptosystem
   - Select images from images/image and Start running cryptosystem.ipynb
   - At the end you will find all results saved in your folders

## Contributing

Contributions are welcome! Feel free to fork the repo, open issues, and submit pull requests.

## Acknowledgements

Special thanks to researchers in cryptography whose work inspired this research.
