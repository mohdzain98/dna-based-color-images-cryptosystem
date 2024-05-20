<a name="br1"></a> 

**Design of DNA Based Color Image Cryptosystem and its Security**

**Analysis.**

ARTICLE INFO

ABSTRACT

The goal of image cryptosystems is to protect image transmission when

there are network adversaries present. To ensure secrecy, images are

subject to encryption to produce unintelligible cipher images; the

techniques used for this process differ significantly from those applied to

text data. The majority of the cryptosystem consider complicated or

confusion–diffusion architectures that changes and permute the values of

the pixels. These frequently entail binary operations like bitwise XOR, plus–

minus, DNA operations, etc; are carried out utilising chaotic maps, each of

which have certain limitations. This paper employs a binary function that

can be applied to both traditional and DNA techniques for coloured natural

images, and can be applied to any kind of image cryptosystem. Each of the

color component of the image follows some steps which start with inter

channel mixing and mix rows which takes the value from initial vector and

mixes it, then Arnold cat map algorithm is applied to shuffle the pixels.

Secondly spiral mixing of the pixel is applied in four different directions to

prevent differential attack. Finally, encoding, substitution, and decoding

based on DNA is carried out. A multiple collapse chaotic map is used to

derive initialization vector, rule charts and DNA substitution map which are

used in encoding and decoding processes. Analysis and experimental

findings shows that our cryptosystem has significant performance and

different metrics shows it can withstand different type of attacks.

*Keywords:*

Color Images Encryption

Color Images Decryption

DNA based Encoding

DNA based Decoding

Experimental Analysis

**1. Introduction**

The ease with which technology is now available, that DNA can be used to store and process

such as smartphones and the Internet, has made information in a way that is similar to how a

sharing images quite simple. These days, a growing computer uses binary digits (bits) to represent data.

number of multimedia images are routinely sent Instead of representing images using traditional pixel

across public networks, and there is an increasing values the image data is converted into a sequence

focus on how to safeguard them. One of the best of DNA nucleotides(A for Adenine, T for Thymine, C

techniques is picture encryption, which converts the for Cytosine and G for Guanine). This process of

important images into noise-like ones that can be converting color images’s pixel values into a DNA

sufficient to protect the hidden images. In contrast to sequence involve mapping the pixel value of a image

text data, images are stored in several formats based to specific sequence of DNA bases. For each color

on their intended use and have highly associated channel i.e. RGB of pixel might be represented by a

properties. Both lossless and lossy formats are used certain DNA bases or combination of bases. Overview

to store natural images. Cryptosystems intended for of Conversion- The value of pixel 0-255 can be

images are distinct from those used for other converted into DNA sequence. Using a set of

purposes since they are made to change pixel predefined rules we can map every two bit of a 8 bit

intensity values and guarantee the removal of binary sequence to one of the DNA nucleotides.

statistical and visual information.

Most image cryptosystems produce quasirandom

*DNA based color images* refer to the encoding of the sequences, wherein the values are predetermined by

color images into DNA sequences called as DNA the initial values supplied to them; these sequences

computing. The basic idea behind DNA computing is serve as the algorithm's symmetric key, by using

1



<a name="br2"></a> 

specific chaotic maps such as the 2D-Infinite Collapse *This research* will use a binary operator (⊗) proposed

Map [1], Logistic Maps [2] 3D chen chaotic sysytem in paper[9]. The operator satisfies the requirements

[3], Hyper Chaotic Systems [4], Logistic Tent System to be utilised with DNA nucleotides and image with

[5], 4D memristive Hyper Chaos [6], Fractional Order pixel depths of 8 and 16 bits.

Hyperchaotic System [7], Fractional Fourier

A new cryptosystem that can be used to safely and

losslessly encrypt color images has been developed

in order to demonstrate the operator's viability with

color images as well. The cryptosystem performs

fairly because it functions well for all color images

with different dimensions and it can withstand a

variety of attacks.

Transform [8] etc.

Most of the cryptosystem follows complex

architectures, Xiuli Chai et al.[3] devised a method of

encryption that creates stochastic sequences using

the Chen chaotic system, which are then used to

form arrays for key stream creation and image

permutation. The two stages of their image

encryption approach are diffusion and permutation.

Mohammad Hossein Moattar and Abolfazl Yaghouti

Niyat [4] proposed a cryptosystem in which color

image is firstly decomposed into R(Red), G(Green) &

B(Blue) channels, permutation is used to mix up all

three components. Using DNA encoding rules

permuted components are converted into DNA

The operator will be applied on RGB colour images in

this study, and its output will be compared to the

output of other cryptosystems that are utilised as

references.

**2. Proposed Work**

matrices. After this, Diffusion is applied. To enhance An encryption and decryption technique utilising an

security they applied second confusion scheme. In operator is presented in the proposed work.

⊗

paper [6] Zhentao Liu et al. suggested use of 4D Measurements of the cryptosystem's performance

memristive hyper chaos with dynamic DNA for colour have also been made using the experimental results

image encryption in which the process to obtain the and analysis. Lastly, a comparison with other

R, G, and B DNA matrices, dynamic encoding is operators, such as XOR and addition-subtraction, is

applied independently to each of the image's three carried out.

channels. For all three matrices, dynamic dispersion

·

The new methods used in the algorithm,

such as inter channel mixing, mix rows, spiral

mixing are straightforward but have the

power to significantly reduce the visual

information in the original images.

The cryptosystem utilizes the modified ACM

algorithm proposed in paper [9], only images

with similar dimensions can be shuffled using

traditional ACM. Unevenly sized

images shuffling is also possible with the

modified ACM.

To improve security and give the encryption

process non-linearity, the suggested

cryptosystem makes use of DNA encoding,

substitution based on DNA with the

proposed operator and DNA decoding.

The cryptosystem has been tested using a

number of statistical and security measures,

together with experimentation and analysis

to improve its performance, to ensure that it

is resistant to a variety of threats.

The comparative analysis shows that

suggested operator outperforms other

operators, such as XOR or (+,−) with the color

images as well, and hence overcomes their

shortcomings.

and dispersal are utilised. Finally, DNA decoding and

component integration resulted in a single encrypted

image. In the paper[8] M. Amine Ben Farah et al.

introduced a cryptosystem with the operations like

DNA encoding, XOR with encoded image and SHA-2

algorithm, substitution, DNA encoding, Fractional

fourier transform applied on encoded image, XOR is

applied on transformed image with Chaotic sequence

obtained from Lorenz system, again substitution,

again Fractional Fourier tranform of this substituted

image, again XOR with chaotic sequence and

obtained substituted image, finally Partial Fourier

Transform and XOR with third chaotic sequence and

ciphered image is obtained.

·

·

·

·

Different authors have used number of algortihms

for color image’s encryption. To improve the

execution of the cryptosystem, they entails several

operations including DNA based encoding, DNA-XOR,

DNA complement, DNA decoding, etc. The majority

of cryptosystems rely on addition–subtraction (+,−)

operations [3] or XOR operations [2,3,4], each of

which has certain drawbacks. The disadvantages of

using XOR in these operation are proved in paper [9].

In the same paper to overcome these limitation a

new binary operator is proposed which is applicable

in both conventional and DNA based cryptosystems.

2



<a name="br3"></a> 

**2.1 Encryption Algorithm Flow Chart**

MCCM chaotic map

240 bit key

0110101010…

Key Scheduling

*RM2*

Initial Vector

Rule Map

DNA Substitution Map

*RM1*

*RM4*

*RM3, DSM*

R

G

B

Inter

Channel

Rows Mix

DR

DG

DB

Rectangle

ACM

Spiral

Mixing

DNA

Encoding

DNA

Decoding

Substitution

Color

Image

cipher

Image

Block Diagram of Ecryption process

**2.2 The Proposed Crptosystem**

using proposed operator with the substitution map.

Lastly, one of the eight DNA decoding guidelines

The cryptosystem makes use of six 2D-Multiple identified by a different rule chart is applied to

Collapse Chaotic Maps, which are produced by the decode the substituted DNA sequences.

Key Scheduling Algorithm from a 240-bit symmetric

*2.2.1 Private Key and Chaotic Map*

key. The nonlinear maps produce quasirandom

sequences in the interval [−1, 1], from which certain

functions are used to construct values and maps such

as DNA substitution maps, rule maps, and initial

vector. The different procedures required to encrypt

the image make use of these variables and mappings.

The suggested approach to create the initial

conditions of the 2D MCCM map, which has been

thoroughly covered in the paper[11], uses

symmetric key with a length of 240 bits.

a

These equations are used to construct the 2D MCCM

as specified in [11]:

The first stage involves intermixing of channels of the

image then utilising the proposed operator ⊗ to

operate the rows of the intermixed images with the

initial array. Using a modified version of Arnold's Cat

Map [10], which can be applied to both square and

rectangular images, the pixels of the composite

image are jumbled. The suggested operator is then

used to further mix the pixels in the directions of

forward spiral row, forward spiral column, reverse

spiral row, and reverse spiral column. These keyless

procedures are necessary to eliminate the link

between pixels and the image's visual perception.

Lastly, a non-linear layer has been added by

combining the suggested operator with the idea of

DNA cryptography. One of the eight encoding

principle identified by the rule chart is applied to the

mixed image that was obtained in the preceding

where the values of the tuning parameters, 푎 and 푏,

are taken from the key. The starting conditions of the

kth dynamical map *M(k)* are represented by (*x0*, *y0*).

The M(1) map is constructed using key and other

maps are constructed using last values of previous

maps.

The operator ⊗ in paper[8] is given as

c ⊗ d = {(c + 1) × (d + 1)mod p} − 1

phases, converting it into a sequence of DNA a) *Initial Vector*

nucleotides. At this stage DNA encoded image is

obtained and these encodings are further substituted With each row of the given image, the initial array

(퐼푉) is utilised to generate a quasirandom sequence

3



<a name="br4"></a> 

that may be operated with ⊗. For an image of size *m* a) Split the R, G, B pixel matrices from the color

*\* n* and bit depth P, the IV-dimension is *1 × m*. The image into three different matrix say red, green,

values are obtained in the following way from 푀(1).

blue. We need to apply all steps to every matrix

inorder to encrypt the image in best way.

InV (1, i) = |xi | × |yi | × (2<sup>ꢂꢃ</sup> − 1)mod (2<sup>ꢄ</sup>),

푥푖 , 푦푖 ∈ 푀(1), 푖 = 1, 2, … 푚

b) *Inter Channel Mixing and Mix Rows*

In this stage, the channels of the image are

intermixied with each other using some equations.

The values from the initial vector (퐼푉), which has

already been generated, are mixed and distributed

throughout the inter mixed image's rows. The IV will

have a dimension of 1 \* 푚 for a plain image 퐼 with

dimensions 푚 \* 푛. To create image 퐼′, the first value

in each row of the matrix is operated with ⊗ using

the appropriate values of the initial array.

b) *Rule Map*

To construct one of the eight DNA encoding-

decoding principle, as listed in Table 1, the rule

mappings are utilised to create a quasirandom

sequence of integers from 1 - 8. Because every two

bits are mapped to one nucleotide, an image of size

푚 \* 푛 of bit depth 푝 will have a rule map of

dimension 푚 × (푛 ⋅ (푝∗2)). The dynamical map is

sclaed to the dimension of image in order to produce The channel are Inter Mixed using these equations:

the rule map. Modular operation is used to

transform the precept, which are preferably selected 푅<sup>ꢅ = 푅</sup>

⊗

<sup>(</sup>G ⊗ B)퐺<sup>ꢅ</sup> = 퐺 <sup>ꢆ</sup>R<sup>′</sup> ⊗ Bꢇ퐵<sup>ꢅ</sup> = 퐵 <sup>ꢆ</sup>R<sup>′</sup> ⊗ G<sup>′</sup>ꢇ

⊗

⊗

from each dimension, to values between 1 and 8.

The value of the first element in each row is

propagated along its neighbor pixels using the

following equation to skew the value in the

neighbouring pixels:

Table 1

DNA Encoding rules

Binary

digits

R1 R2 R3 R4 R5 R6 R7 R8

퐼′ (g, h) = I(g, h) ⊗ I′ (g, h − 1);

ℎ = 2, 3 … , 푛; ∀푔 ∈ {1, … , 푚}

00

A

C

G

T

A

G

C

T

C

A

T

G

A

T

C

T

A

G

G

T

A

C

T

C

G

A

T

G

C

A

01

10

11

The process of inter channel mixrows will be applied

on the image received after splitting color image into

matrices i.e. R, G, B, respectively as.

G

C

C) *DNA Substitution Map*

푅′ = 푖푛푡푒푟푀푖푥(푟, 푔, 푏)

퐺′ = 푖푛푡푒푟푀푖푥(푔, 푅<sup>ꢅ</sup>, 푏)

퐵′ = 푖푛푡푒푟푀푖푥(푏, 푅<sup>ꢅ</sup>, 퐺′)

Using the proposed operator ⊗, the DNA

Substitution chart—a pseudorandom string of DNA

nucleotides (퐴, 퐶, T, 퐺)—is utilised to work with the

mixed image which is DNA encoded. With size of

*m\*n* an image with a *p* bits pixel depth will be

enlarged to a DNA sequence of size 푚 × (푛.(푝∕2))

which the DNA Substitution Chart's size ought to

correspond to. The instruction for constructing a

DNA Substitution chart is described in paper[9].

푚푖푥푒푑퐼푚푎푔푒푅푒푑 = 푚푖푥푅표푤푠(푅<sup>ꢅ</sup>, 퐼푉

푚푖푥푒푑퐼푚푎푔푒퐺푟푒푒푛 = 푚푖푥푅표푤푠(퐺<sup>ꢅ</sup>, 퐼푉

푚푖푥푒푑퐼푚푎푔푒퐵푙푢푒 = 푚푖푥푅표푤푠(퐵<sup>ꢅ</sup>, 퐼푉)

)

)

After obtaining the necessary maps, the images can

be encrypted. The decryption procedure also

requires the same maps. The next subsection

provides a description of the encryption procedure.

2\.2.2 *Encryption*

An image in color is changed into an

incomprehensible distorted image during the inscribe

process. The following steps make up the encryption

process:

The algorithm for mixrows is same as defined in

paper[9]

c)*Rectangular ACM*

4



<a name="br5"></a> 

d) *Spiral Mixing*

Arnold’s Cat Map(ACM)[10] is a disorganised chart

that is used to rearrange the given image. When a 2D

array has dimensions of N\*N, ACM is expressed as:

In order to prevent differential attacks on the image,

this stage job is to propagate the pixel’s value

information throughout it in spiral manner. This

procedure helps remove the image's visual

information even if it doesn't require an encryption

key. Let 푚 and 푛 represent the input image 퐼's row

and column counts, respectively. There are four parts

to the mixing step.

푔′

ℎ′

2

1

1

1

푔

ℎ

= ꢈ

ꢉ

(푚표푑 푁)

where g, h, g′, h′ ∈ [1, 푁]. (g, h) and (g′, h′)

represents, respectively, the original and scrambled

images' position indices. The drawback of employing

directly ACM algorithm is that its definition is limited

to identical size images. Nonetheless, uneven

proportions are a regular occurrence in both natural

and medical imaging. The image is separated into

minimal squares, and ACM is applied independently

to each of them in order to fix this problem.

·

Forward Spiral Row Mixing

The mixing is started from I[0,0] and goes to

[m,n] in spiral fashion. Every element value is

propagated

Let m\*n represent an image's dimension. To divide it

into as few squares as feasible, the length of each

square's side will be 푁 = min(푚, 푛). Therefore, the

proportion of squares that are obtained will be 훼 =

⌈max(푚, 푛)/푁⌉. If m and n are distinct combination

of one another, then 퐿 = 푁 − (max(푚, 푛) mod 푁)), an

additional length L pixel remains with the longer side.

There will be some overlap between the squares in

order to accommodate this additional space within

the image. There can be a maximum of 훼 − 1

overlaps to ensure almost equal overlapping; the

duration of the last overlap i.e. eta=⌊퐿∕(훼 − 1)⌋. The

last overlap modifies the excess portion that is still

present.

Given that 푚 < 푛 for an image, 푁 = 푚 and the initial

value of I[1,1] of first square, a = 1, b = 1. For the

initial 훼 − 1 squares, the subsequent values will be

ascertained as b′ = b + 푁 − eta. To accommodate the

extra square that remains, the 푦-coordinate for the

last square will be 푦−푁. For images when 푚>푛, same

computations can be performed by transposing the

푚 and 푛 values. Next, ACM is applied to each

separated square individually. On the other hand,

ACM is regular; that is, the original image is restored

following a definite number of iterations. The image

푁's dimension determines this regularity. To achieve

a well-shuffled image and avoid regularity, ACM will

be applied to the image itself t number of times,

where t is the highest prime number less than or

equal to (푁/2), 푁 being the size of each image.

·

Forward Column Mixing

As forward row, in this step we move in

column from M[0][0] to M[0][n]. The values

of each element is propagated in the

columns.

The process of applying ACM on images will be

applied on the image received after mixrows

operation on every matrix i.e. R, G, B, respectively as

푖푚푎푔푒퐴퐶푀푅 = 푎푝푝푙푦퐴퐶푀(푚푖푥푒푑퐼푚푎푔푒푅푒푑)

푖푚푎푔푒퐴퐶푀퐺 = 푎푝푝푙푦퐴퐶푀(푚푖푥푒푑퐼푚푎푔푒퐺푟푒푒푛)

푖푚푎푔푒퐴퐶푀퐵 = 푎푝푝푙푦퐴퐶푀(푚푖푥푒푑퐼푚푎푔푒퐵푙푢푒)

5



<a name="br6"></a> 

·

Backward Spiral Row Mixing

Reverse spiral rows mixing is the reverse of

rows mixing in which we start from M[m,n]

and go towards M[0,n]

The values are operated with each other.

[a]

[b]

[c]

[d]

**Figure:** a – d are the spiral mixing operations described

ahead.

e) *DNA Encoding*

One of the encoding principle transforms each pixel

value in the image into

a

string of DNA

nucleotides. At a time, two bits are extracted and

transformed into a nucleotide by applying the rule

found in the obtained rule map. As a result, a string

half of length i.e. (p/2) nucleotides per pixel will be

created from an image with a bit depth of 푝-bit.

The encoding process of an image will be applied on

the image received after mixing on every matrix i.e.

R, G, B, respectively as

·

Backward Spiral Column Mixing

In this step we go backward from M[m,n]

푒푛푐표푑푒푑푅 = 푒푛푐표푑푖푛푔(푖푚푔푚푖푥푅)

towards zero in column way

푒푛푐표푑푒푑퐺 = 푒푛푐표푑푖푛푔(푖푚푔푚푖푥퐺

)

The element value is propagated from the

푒푛푐표푑푒푑퐵 = 푒푛푐표푑푖푛푔(푖푚푔푚푖푥퐵)

very last value to the first value of the matrix

25

12

76

54

78

98

65

68

31

CACA CGCG ACAC

CCTT ACTG ACCT

ACTT ATAG CATA

**Figure:** Illustration of Encoding process

f) *Substitution*

The DNA sequences 퐷 that are received after the

image has been encoded are substituted in the

replacement process with ⊗ using the sequence that

is obtained from the DNA replacement Map (DSM).

The mapping of 훴 = {퐴, 퐶, 푇, 퐺} to S = {0, 1, 2, 3} for

DNA operation requires the selection of one among

eight encoding rules from the values derived from

the rule map (RM). Eight options will appear when

⊗ is used, one for each of the encoding rules that

are applied to a nucleotide’s pair; yet, some pairs

may give identical outcomes for distinct rules. Since

of this, the dynamical DNA sequence that is left over

after this procedure is resistant to algebraic attack

since it depends on the rule’s table as well as the

elements of the DNA substitution table, both of

which are necessary to decrypt the cryptic image

and determine the original key.

These four steps are required to be processed on the

image one after other in continuation on every image

taken into consideration. The process of mixing

image will be applied on the image received after

applying ACM on every matrix i.e. R, G, B,

respectively as

푖푚푔푚푖푥푅 = 푠푝푖푟푎푙푀푖푥푖푛푔(푖푚푎푔푒퐴퐶푀푅)

푖푚푔푚푖푥퐺 = 푠푝푖푟푎푙푀푖푥푖푛푔(푖푚푎푔푒퐴퐶푀퐺)

푖푚푔푚푖푥퐵 = 푠푝푖푟푎푙푀푖푥푖푛푔(푖푚푎푔푒퐴퐶푀퐵)

6



<a name="br7"></a> 

The following represents the outcome of DNA operations that correspond to the various encoding rules:

A

A

C

G

T

C

C

T

G

G

A

T

T

A

A

C

G

T

C

C

T

G

G

A

T

T

A

G

A

T

C

A

C

G

T

G

T

T

A

C

T

C

T

G

A

V

G

T

T

A

C

G

T

T

A

C

G

T

T

A

C

G

T

C

T

A

C

G

T

G

A

G

C

A

G

C

A

G

C

A

G

C

A

A

G

A

G

A

G

A

G

T

C

C

C

C

R1

R2

R3

R4

A

G

A

T

C

A

C

G

T

G

T

T

C

T

A

G

A

T

C

A

C

G

T

G

T

T

C

T

A

T

C

G

A

T

G

C

T

T

A

C

G

T

A

T

C

G

A

T

G

C

T

T

A

C

G

T

A

C

G

T

A

C

G

T

A

C

G

T

A

C

G

T

G

C

A

G

C

A

G

C

A

G

C

A

A

G

A

G

A

G

A

G

C

C

C

C

R5

R6

R7

R8

Both the rule map (RM) and the DNA substitution process is the final cipher image. The algorithm for

map (DSM) must be used during the substitution DNA decoding to decimal value is given in paper[9].

procedure. Following encoding on each matrix, R, G,

The process of decoding image will be applied on the

image received after applying substitution on every

matrix i.e. R, G, B, respectively as

and B, respectively, the process will be applied to the

image received as

푠푢푏푠푡푖푡푢푡푒푑푅 = 푠푢푏푠푡푖푡푢푡푖표푛(푒푛푐표푑푒푑푅, 퐷푆푀, 푅푀)

푠푢푏푠푡푖푡푢푡푒푑퐺 = 푠푢푏푠푡푖푡푢푡푖표푛(푒푛푐표푑푒푑퐺, 퐷푆푀, 푅푀)

푠푢푏푠푡푖푡푢푡푒푑퐵 = 푠푢푏푠푡푖푡푢푡푖표푛(푒푛푐표푑푒푑퐵, 퐷푆푀, 푅푀)

푓푖푛푎푙퐶푅 = 푑푒푐표푑푖푛푔(푠푢푏푠푡푖푡푢푡푒푑푅)

푓푖푛푎푙퐶퐺 = 푑푒푐표푑푖푛푔(푠푢푏푠푡푖푡푢푡푒푑퐺)

푓푖푛푎푙퐶퐵 = 푑푒푐표푑푖푛푔(푠푢푏푠푡푖푡푢푡푒푑퐵)

Encoding, substitution and decoding all three

processes uses three different rule maps. As shown

in flow chart substitutionuses rule map 3.

Merging all three decoded matrix into final cipher

image. The function merge is inbuilt function of

python’s opencv library.

g) *Decoding*

푐푖푝ℎ푒푟퐼푚푎푔푒 = 푚푒푟푔푒(푓푖푛푎푙퐶푅, 푓푖푛푎푙퐶퐺, 푓푖푛푎푙퐶퐵)

One of the eight decoding rules is used to the chaotic

DNA sequences in order to transform them to The cipher image that is produced is noise-like and

decimal values; the rule map's sequences determine random, meeting the criteria for a secure image. This

which rule is applied. Then, a single cipher color cipher picture is resistant to cryptographic attacks,

image is created by combining the three matrices— has low correlation, and high entropy. The

cipher red, cipher green, and cipher blue obtained decryption procedure and security analysis are

after decoding. The outcome of this decoding displayed up front.

**2.3 Decryption Algorithm Flow chart**

MCCM chaotic map

240 bit key

1101010010…

Key Scheduling

DNA Substitution Map

*RM2* Rule Map

Initial Vector

*RM4*

*RM3, DSM*

*RM1*

DR DNA

DG Encoding

DB

Reverse

DNA

Substitution

DNA

Decoding

Reverse

Spiral

Mixing

Reverse

Rectangle

ACM

Reverse

Inter

Channel

Rows Mix

R

G

B

Decrypted

Image

Cipher

Image



<a name="br8"></a> 

*2.3.1 Decryption*

apply reverse modified ACM followed by reverse mix

rows and reverse inter channel mixing. At last we

After obtaining cipher image we proceed for its again merge RGB into final decrypted image. The

decryption as on reciever side it needs to be main point of decryption is to use the same key as it

decrypted. The encryption and decryption processes was used in encryption for getting 2D MCCM, initial

are inverted. The steps that are taken in last stages vector, rules maps and DNA substitution map. The

will now be taken in initial stages as we move in

reverse fashion. As presented in flow chart, we again

split the cipher images into R, G, B components and

applies DNA encoding. If we compare encyption and

decryption flow chart we can see we are dealing with

DNAs in the later stages in encryption but in

decryption we are using DNAs in initial stages. At first

we do DNA encoding then we apply reverse DNA

substitution using reverse substitution rules followed

by DNA decoding. Encoding and Decoding have the

same algorithm but substitution have minor change,

it uses reverse substitution rules. After decoding we

go for reverse spiral mixing which internally follows

reverse fashion of mixing in encryption. Then we

steps are shown ahead.

a) Splitting final crypted image into R for Red, G for

Green and B for Blue components.

b) *DNA encoding*

As a first step now DNA encoding is performed on

the R, G, B, matrices received after splitting. There is

no change in the algorithm of encoding.

c) *Reverse Substitution*

For reverse substitution we define reverse encoding

rules between DNA nucleotides.

Reverse Substitution Rules

A

A

C

G

T

C

G

A

T

G

C

T

A

G

T

T

G

C

A

A

C

G

A

T

C

G

C

T

A

G

T

T

G

C

A

A

C

T

A

G

C

A

C

G

T

G

T

G

C

A

T

A

G

A

G

C

C

T

G

C

A

G

A

V

G

T

T

C

T

A

G

A

C

G

T

A

C

G

T

A

C

G

T

A

C

G

T

G

A

T

A

C

G

T

C

C

R1

R5

R2

R6

R3

R7

R4

R8

A

C

T

A

G

C

A

C

G

T

G

T

G

C

A

T

G

A

T

A

G

A

T

C

T

G

C

A

G

A

C

G

T

T

C

T

A

G

A

A

G

C

A

C

G

T

A

G

G

C

A

T

T

T

C

G

T

A

A

C

C

A

C

G

A

A

G

G

C

T

T

C

T

T

G

G

T

A

C

G

T

A

C

G

T

A

C

G

T

A

C

G

T

C

C

C

d) *DNA Decoding*

The decoding process has the same algorithm as we

have used in encryption process defined in paper[9].

For further steps of decryption we reverse the

operator ⊗ in decryption algorithm as

(a ⊗<sup>R</sup> b) ꢆ(푎 + 1) ∗ 푚표푑푖푛푣(푏 + 1, 푝)ꢇ푚표푑 푝 – 1

\=

e) *Reverse Spiral Mixing*

In reverse mixing we apply mixing in reverse fashion.

Rev mixing have same four steps but in reverse

order.

The algorithm will require to use DSM map and rule

map generated at the time of decryption with the

same key and same chaotic maps.

·

Backward Spiral column mixing

In the process of decryption we move from

last operation towards first.

푠푟푒푑 = 푟푒푣푆푢푏푠푡푖푡푢푡푖표푛(푒푟푒푑, 퐷푆푀, 푅푀)

푠푔푟푒푒푛 = 푟푒푣푆푢푏푠푡푖푡푢푡푖표푛(푒푔푟푒푒푛, 퐷푆푀, 푅푀)

푠푏푙푢푒 = 푟푒푣푆푢푏푠푡푖푡푢푡푖표푛(푒푏푙푢푒, 퐷푆푀, 푅푀)

8



<a name="br9"></a> 

f) *Reverse Modified ACM*

In reverse modified ACM we just reverse the loops of

ACM algorithm.

g) *Reverse Inter channel mixing and Mix Rows*

The rows of the images obtained after applying

reverse modified ACM are treated with the same IV

and reverse of operator ⊗.

·

Backward Spiral row mixing

After applying mix rows, the channel are reversly

diversed to get the original r, g, b matrices. We

defined inter mix equations in such a way that they

can be reversed easily. Below equations shows the

reverse procedure. For reverse intermix we start with

getting B channel.

·

Reverse Forward Spiral column mixing

퐵 = 퐵<sup>ꢅ</sup> ⊗<sup>R</sup> ꢆR ⊗ G ꢇ

′

′

퐺 = 퐺<sup>ꢅ</sup> ⊗<sup>R</sup> ꢆ푅<sup>′</sup> ⊗ Bꢇ

푅 = 푅<sup>ꢅ</sup> ⊗<sup>R</sup> ( 퐺 ⊗ B)

·

Reverse Forward Spiral row mixing

h) Merge the decrypted R, G, B matrices to get final

decrypted image.

**2.4 The Final Results**

a)

9



<a name="br10"></a> 

(1)

(5)

(2)

(6)

(3)

(4)

b)

c)

d)

(7)

(8)

**Fig 2.** The results of sensitivity of a key of (1) Lenna, (2) Baboon,

(3) pepper and (4) Lake. (5) to (8) are decrypted images of

respective color images.

*3.2 Histogram Analysis*

A histogram describes how the image's pixel intensity

values are distributed. The histograms of the plain

and encrypted images are present in Figure 3. It is

evident that the cypher image's pixel values are

consistently distributed throughout the interval [0,

255], in stark contrast to the plain image's pixel

values which follows uneven trend. Single image

shows the histogram for all three channels.

Additionally, we use histogram variances to

statistically test the regularity of histograms.

Mathematically, histogram variances are defined as :

**Figure 1.** Original, en-crypted & de-crypted image of a) Lenna

b) baboon c) Pepper d) Lake

**3. Experimental Results and Analysis**

Numerous colour images were encrypted, and the

cipher image’s characteristics were examined in

order to evaluate the cryptosystem's performance.

푚

푛

∑

∑

(I(x, y) − μ)<sup>ꢋ</sup>

*3.1 Space and Sensitivity of key*

푥 = 1 푦 = 1

푚 ∗ 푛

variance(I) =

The encryption's resistance to brute-force attacks is

tested by the length of a key. In aspect of

cryptography, a high degree of security requires a

key length that is greater than 2<sup>100</sup>. In our proposed

cryptosystem, The key used has a secured key space

of 2<sup>240</sup> and is 240 bits long. The key used in

encryption and decryption is same because it is

symmetric.

Colour images have low variance since they only

have eight bits per pixel and have a non-uniform

distribution of pixels in each R, G, and B component

image. Table 5 displays natural images with variances

less than 3000, while the matching cipher images of

the corresponding Re, Gr, and Bl components have

variances greater than 5400.

High sensitivity to the key is a fundamental

requirement for the encryption system; otherwise, it

will not be able to obtain the right decryption image,

even in the event that there is a very slight shift in

the input key. A small modification to the initial key

will result in a unseen new decryption image for a 1.

perfect encryption scheme. In our experiment, we

each alter four RGB colour images by merely varying

the lone bit of the initial key. The decryption

outcomes are displayed, along with the notable

variations between the original (256 x 256) and

decrypted images.[3]

2\.

10



<a name="br11"></a> 

Better pixel intensities are distributed when 휒2 is

smaller. From the 휒2-Inverse Cumulative Distribution

Function [12], we can determine the crucial values

(upper limit) of 휒2 for 8-bit pictures with level of

significances 훼 = 20. These values are 273.79. Since

the produced cipher images, 휒2 values are smaller

than the crucial values, The hypothesis holds for

coming in the range of 20% significance, indicating a

uniform variation of pixel magnitude across all of

them.

3\.

4\.

*3.3 Entropy Analysis*

The uncertainty measurement is provided by

entropy. The entropy ɛ(푆) is calculated for a source 푆

that contains 푛 symbols, each with a chance of

presence of 푝푖 as[16]

**Figure 3. 1**) plain and cipher histogram of Lenna **2**) plain and

cipher histogram of baboon **3**) plain and cipher histogram of

pepper **4**) plain and cipher histogram of lake.

Histogram of color images are having different

graphs and bins of different sizes of all three color

channels but cipher image have all three channel

equally distributed on the graph.

푛

푖 = 1

1

ɛ(S) = ꢗ

pi ∗ log (

)

pi

Entropy of plain images is low, but that of cipher

images should be close to 7 high bits/symbol. Table

ahead demonstrates that the crypted images

produced by the suggested technique have an

entropy larger than 7.9950 , indicating a high degree

of randomness.

*3.3 Chi-Square Analysis*

The encrypted image's histogram offers a visual

depiction of the data, while the chi-Square test offers

a statistical analysis of the variation of pixel intensity

readings. The chi<sup>2</sup> readings for an image of size 푚 × 푛

and pixel depth 푝 bits can be computed as [15]

*3.4 Correlation Analysis*

<sub>2ꢌ − 1 ꢆꢍꢎꢏ – ꢍꢐꢇ</sub>ꢑ

퐿 = 0 ꢒ<sub>ꢓ</sub>

Pixels have a strong correlation in plain images,

however in crypted images, this correlation should

be reduced. One among key statistical features of

images is the correlation between neighbouring

pixels, which indicates the limit of correlation

between pixel values in adjacent positions in the

image. The greater the efficiency of the developed

encryption system, there is low correlation between

the neighbouring pixels of the encrypted picture.

Next, 10000 pixels are chosen from both the

encrypted and original images in order to examine

χ<sup>2</sup> =

∑

<sup>where xꢔ</sup><sub>ꢕ</sub> is the frequency that has been observed

for the 퐿th intensity, where 퐿 is between 0 2<sup>ꢌ</sup>− 1,

<sup>and 푥</sup>ꢖ <sup>is the frequency expected, which has been</sup>

computed as:

ꢌ

푥 = 2 − 1/(푚 ∗ 푛)

ꢖ

The table below displays the values for 휒2 for each of

the cipher images.

Table 4: chi-Square of cipher images

the pixel correlation. The coefficient of correlation

휌

is calculated as [17]

Image

Chi-Square

Cipher Image

Red Green

ꢚ

(∑ (x − xꢙ)(y − yꢙ))

ꢘꢛꢃ

ꢘ

ꢘ

Blue

휌(x, y) =

ꢜ∑ <sup>(x</sup>ꢘ <sup>− xꢙ)ꢋ ꢜ</sup>

ꢚ

ꢘꢛꢃ

∑<sup>ꢚ</sup> (<sup>y</sup><sub>ꢘ</sub> <sup>− yꢙ)ꢋ</sup>

ꢘꢛꢃ

Lenna

Baboon

Pepper

Lake

267\.47 230.02

232\.32 268.43

246\.94 216.14

257\.71 235.82

270\.61

267\.44

246\.22

272\.62

푛 푥

Consider two vectors of size , and , respectively,

with means̄and̄. -1 ≤ ≤ 1. The original images

have a high positive correlation when the values are

close to 1, while the cipher images benefit from a low

correlation when the values are close to 0.

푦

푥

푦

휇

11



<a name="br12"></a> 

Table 5: Variances of plain and cipher images

Image

Plain Image

Cipher Image

Red

Red

Green

Blue

Green

Blue

Lenna

Baboon

Pepper

Lake

2412\.122 2757.27 1132.811

2889\.991 1972.825 3383.311

1995\.873 2556.325 1881.152

1954\.858 2572.911 2548.812

5464\.239 5429.818 5451.349

5459\.752 5444.902 5457.723

5469\.795 5444.531 5485.011

5469\.953 5498.777 5421.363

Table 6: Entropy of plain and cipher images

Image

Entropy

Plain Image

Red

Cipher Image

Green

Blue

Red

Green Blue

Lenna

Baboon

Pepper

Lake

7\.268

7\.683

7\.271

7\.331

7\.597

7\.381

7\.527

7\.626

6\.971

7\.682

7\.109

7\.341

7\.9951

7\.9972

7\.9972

7\.9971

7\.9974 7.9969

7\.9973 7.9971

7\.9975 7.9972

7\.9974 7.9973

*3.5 NPCR and UACI*

Table 7: NPCR of cipher images

Images

NPCR

N: Number of P: Pixel C: Chaging R: Rate and U:

Unified A: Averaged C: Changing R: Rate are crucial

metrics that assess a cryptosystem's resilience to

differential attacks. Let X and X′ represent the *c*ipher

images that are generated following the *e*ncryption

of a plain *i*mage of 푚 \* 푛 size and a pixel depth of 푝

i.e. X, and modifying one value of pixel of the initial

image respectively and encrypting again is X’. The

N.P.C.R. and U.A.C.I. are computed using the

following formulas..

NPCR critical value = 99.534077

Red

Green Blue

Lenna

Baboon

Pepper

Lake

99\.6086

99\.6062

99\.5971

99\.5994

99\.6018 99.5924

99\.6062 99.6154

99\.5986 99.6121

99\.5986 99.6185

Table 8: UACI of cipher images

1, 푖푓 푋(푢, 푣) ≠ 푋<sup>ꢅ</sup>(푢, 푣)

0, 푖푓 푋(푢, 푣) = 푋<sup>ꢅ</sup>(푢, 푣)

Images

UACI

퐷(푢, 푣) = ꢝ

UACI range: 33.1593 - 33.7676, α = 0.005

퐷(푢, 푣)

Red

Green Blue

푁

푁

= ꢗ

∗ 100%

Lenna

Baboon

Pepper

Lake

푚 ∗ 푛

33\.4419

33\.3902

33\.4408

33\.4693

33\.3474 33.5046

33\.4853 33.2344

33\.4035 33.4181

33\.3303 33.4835

ꢞ,ꢌ

|푋(푢, 푣) − 푋<sup>ꢅ</sup>(푢, 푣)|

2<sup>ꢌ</sup> − 1 · (푚 ∗ 푛)

푈

= ꢗ

ꢞ,ꢌ

where u = 1, 2, …, 푚; v = 1, 2, …, 푛 & |t| denotes the

absolute value of t.

*3.6 **P**eak **S**ignal to **N**oise **R**atio*

In [13], the optimal settings for UACI and NPCR have **P**eak **S**ignal-to-**N**oise **R**atio (P.S.N.R) is a useful metric

been examined. The greatest value of Fmax for a for evaluating the effectiveness of a cryptosystem

picture with dimensions of m \* n and bit depth of p against noise attacks and occlusion. Let Pl be the

will be F = 2<sup>p</sup> – 1. The table ahead shows N.P.C.R and plain image that is encrypted to create the crypted

U.A.C.I values observed for different images

image Z = En(Pl, key) using a key and an encryption

12



<a name="br13"></a> 

technique En. The crypted image is subjected to

noise attack or occlusion to produce the image Z′.

This image is then decrypted using the matching

decoding algorithm Dc and the same key to provide

the decrypted *i*mage D = Dc(Z′, key). The PSNR is

defined in Equation as follows[14] if the image has

dimensions of m \* n and bit depth of p.

Table 9: PSNR of cipher images

Images

PSNR

Red

Green Blue

Lenna

Baboon

Pepper

Lake

8\.27

8\.76

8\.84

8\.81

8\.25

8\.76

8\.76

8\.81

8\.22

8\.72

8\.76

8\.8

(2<sup>ꢌ</sup> − 1)<sup>ꢋ</sup>

푃. 푆. 푁. 푅 = 10 ∗ log<sub>ꢃꢟ</sub>

ꢠ

ꢥ

1

∑<sup>ꢣ</sup> ∑<sup>ꢡ</sup><sub>ꢢꢛꢃ</sub>{퐷(푢, 푣) − 휌(푢, 푣)}<sup>ꢋ</sup>

**4. Comparative Analysis**

푚 ∗ 푛 ꢤꢛꢃ

Using the same experimental setup and

cryptosystem, the operator ⊗ defined in paper[9]

has been compared with other operators such as

XꚚR and (+/−) on several color images. The NPCR and

An image needs a high PSNR—ideally more than

30dB— to be visually discernible. However, because

of their random noise like properties and high mean

squared error in the associated pixels, cipher images

have low PSNR values—typically less than 10dB.

UACI tests demonstrate that the operator

⊗

outperforms XꚚR in terms of resistance to

differential attacks, since the operator ⊗ yields

satisfactory results whereas XOR fails horribly in

these tests using with the suggested cryptosystem.

While the (+,−) operator yields fairly similar values

for UACI and NPCR but also fails in some

circumstances as well like for one channel fails and

for other two channel give satisfactory result.

[a]

[b]

However, the main disadvantage of (+/−) operator is

removing the pixel correlation and the suggested

operator works significantly better in these

situations. These findings support the assertion that

the operator ⊗ outperforms the standard operators

typically employed in traditional image

cryptosystems. Table 11 alse provides an overview of

various existing cryptosystem’s performances and

compare our values with theirs values.

[c]

[d]

**Figure-4.** Occlusion attack – [a] Original Lenna image, [b]

encrypted image, [c] occlusion attack on encrypted image, and

[d] decrypted image.

Table 10 : Pixel correlation values between plain and cipher images.

Image Direction

Pixel correlation coefficient

Plain Images

Red Green Blue

Cipher Images

Green

Red

Blue

Lenna Horizontal

Vertical

0\.9944 0.9979 0.9845 -0.0002 -0.0036 -0.0014

0\.9956 0.9935 0.9857

0\.9433 0.9268 0.8702

0\.0002

0\.0094

0\.0075

0\.0668

-0.0004

-0.0064

Diagonal

Baboon Horizontal

Vertical

0\.9126 0.8962 0.9531

0\.6595 0.7867 0.8541

0\.0031

0\.0064

0\.0084

-0.0169

-0.0076 -0.0852

Diagonal

0\.7652 0.6806 0.8443 -0.0384

0\.0239

0\.0259

Pepper Horizontal

Vertical

0\.9285 0.9465 0.9486 -0.0042

0\.9725 0.9772 0.9784 -0.0053

0\.0092

0\.0041

-0.0086

0\.0045

Diagonal

0\.9119 0.9455 0.8801 -0.0345 -0.0127 -0.0057

13



<a name="br14"></a> 

Lake

Horizontal

Vertical

Diagonal

0\.9322

0\.9794

0\.8935

0\.9589

0\.9868

0\.9496

0\.9647 0.0042 0.0082 0.0062

0\.9863 -0.0085 -0.0063 0.0058

0\.9458

-0.002 0.0074 -0.0069

Table 11: Comparative Analysis

Image Cryptosystem Encryption Process

Entropy

Gr

NPCR

UACI

Gr

Re

Bl

Re

~

Gr

~

Bl

~

Re

~

Bl

~

Rgb

7\.268 7.597 6.971

color

Original

~

image

rgb color

image

proposed

using ⊗

using XOR

using + -

7\.997 7.997 7.996 99.608 99.601 99.592 33.443 33.347 33.506

7\.997 7.996 7.996 48.416 48.416 48.428 17.39 17.38 17.38

7\.997 7.997 7.996 96.432 99.182 99.424 33.401 33.421 33.405

rgb color

image

Refr [3]

Refr [7]

Dynamic DNA

encryption and

chaos

7\.997 7.996 7.997 99.6 99.6 99.61 33.56 33.46 33.45

7\.997 7.997 7.996 99.6 99.62 99.58 33.47 33.44 33.57

rgb color

image

fractional-order

hyperchaotic

system

and DNA-Coding

rgb color

image

Refr [8]

Refr [4]

fractional fourier

transform and DNA

sequence

7\.997 7.997 7.999 99.56 99.55 99.57 33.41 33.44 33.45

operation

rgb color

Image

hybrid chaotic

system

7\.997 7.997 7.997 99.64 99.6 99.61 33.47 33.36 33.44

and DNA

sequences

The values that fail the crucial test are shown by underlining in the NPCR and UAC

shown in a suggested cryptosystem that includes

DNA encoding, substitution, and decoding after a

unique mixing operation that makes advantage of

the distributive and commutative features of ⊗. The

cryptosystem also demonstrate decryption process

and algortihms where the changes are required for

decryption. In contrast to the existing cryptosystems,

which employ intricate encryption techniques or a

number of rounds for confusion-diffusion, the

proposed cryptosystem just requires basic

computational operations. Analysis and experimental

findings demonstrate that the suggested work

satisfies the specifications needed for a safe images

cryptosystem. Comparative analysis with other

commonly used operators with color images shows

Color images can be successfully encrypted using the

proposed cryptosystem and the operator ⊗ as we

have proposed as well as plain images as shown in

paper [8] . The comparison of the findings shows that

the proposed cryptosystem performs comparably to

the the most modern cryptosystems currently in use.

**5. Conclusion**

In our paper, we design a cryptosystem using a ⊗

operator for color images. The cryptosystem follows

the properties of operator ⊗ which says it forms an

abelian group due to which it is used for encryption

as well as decryption. This operation's application is

14



<a name="br15"></a> 

that operator ⊗ outperforms them regarding color Sci. Technol. J. Sel. Areas Telecommun. (JSAT) 1 (2) (2011) 31–38.

http://refhub.elsevier.com/S0923-5965(21)00173-9/sb26

images encryption.

[14] X. Chai, H. Wu, Z. Gan, D. Han, Y. Zhang, Y. Chen, An efficient

approach for encrypting double color images into a visually

meaningful cipher image using 2D compressive sensing, Inform.

Sci. 556 (2021) 305–340, http://dx.doi.org/10.

1016/j.ins.2020.10.007.

**6. References**

[1] W. Cao, Y. Mao, Y. Zhou, Designing a 2D infinite collapse map

for image encryption, Signal Process. (2020) 107457,

http://dx.doi.org/10.1016/j.sigpro. 2020.107457

[15] X. Wang, L. Feng, H. Zhao, Fast image encryption algorithm

based on parallel computing system, Inform. Sci. 486 (2019) 340–

358, http://dx.doi.org/10.1016/ j.ins.2019.02.049.

[2] Q. Zhang, L. Guo, X. Wei, Image encryption using DNA

addition combining with chaotic maps, Math. Comput. Modelling

52

(11–12)

(2010)

2028–2035,

http://dx.doi.org/10.1016/j.mcm.2010.06.005.

[16] ] X. Chai, H. Wu, Z. Gan, D. Han, Y. Zhang, Y. Chen, An

efficient approach for encrypting double color images into a

visually meaningful cipher image using 2D compressive sensing,

Inform. Sci. 556 (2021) 305–340, http://dx.doi.org/10.

1016/j.ins.2020.10.007.

[3] X. Chai, X. Fu, Z. Gan, Y. Lu, Y. Chen, A color image

cryptosystem based on dynamic DNA encryption and chaos,

Signal Process. 155 (2019) 44–62, https://doi.org/10.

1016/j.sigpro.2018.09.029.

[17] J. Benesty, J. Chen, Y. Huang, I. Cohen, Pearson correlation

coefficient, in: Noise Reduction in Speech Processing, Springer,

2009, pp. 1–4, http://dx.doi.org/10. 1007/978-3-642-00296-0\_5

[4] Abolfazl Yaghouti Niyat & Mohammad Hossein Moattar, Color

image encryption based on hybrid chaotic system and DNA

sequences https://link.springer.com/article/10.1007/s11042-

019-08247-z

[5] X. Chai, J. Bi, Z. Gan, X. Liu, Y. Zhang, Y. Chen, Color image

compression and encryption scheme based on compressive

sensing and double random encryption strategy, Signal Process.

176 (2020) 107684, http://dx.doi.org/10.1016/j.sigpro.

2020\.107684.

[6] Z. Liu, C. Wu, J. Wang, Y. Hu, A color image encryption using

dynamic DNA and 4-D memristive hyper-chaos, IEEE Access 7

(2019)78367–78378,

https://doi.org/10.1109/ACCESS.2019.2922376.

[7] H. Dong, E. Bai, X.-Q. Jiang, Y. Wu, Color image compression-

encryption using fractional-order hyperchaotic system and DNA

coding, IEEE Access

8

(2020) 163524– 163540,

https://doi.org/10.1109/ACCESS.2020.30223 98.

[8] M.A.B. Farah, R. Guesmi, A. Kachouri, M. Samet, A novel

chaos based optical image encryption using fractional Fourier

transform and DNA sequence operation, Opt. Laser Technol. 121

(2020) 105777.

[9] P. Mishra, C. Bhaya, A.K. Pal, A.K. Singh, A novel binary

operator for designing medical and natural image cryptosystems,

Signal Process. Image Commun. 98 (2021) 116377.

https://doi.org/10.1016/j.image.2021.116377

[10] G. Petersen, Arnold cat map survey, Math 45 linear algebra,

http://refhub.elsevier.com/S0923-5965(21)00173-9/sb20

[11] Yang, C., Wei, X. and Wang, C., 2021. S-Box design based on

2D multiple collapse chaotic map and their application in image

encryption. *Entropy*, *23*(10),p.1312.https://www.mdpi.com/1099

-4300/23/10/1312

[12] chi-square inverse cumulative distribution function

\-

MATLAB chi2inv Math Works India, 2021

–

,https://in.mathworks.com/help/stats/chi2inv.html (Accessed 13

May 2021)

[13] ] Y. Wu, J.P. Noonan, S. Agaian, et al., NPCR and UACI

randomness tests for image encryption, Cyber J.: Multidiscip. J.

15

