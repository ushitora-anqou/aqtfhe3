# aqTFHE3

[@nindanaoto](https://github.com/nindanaoto)による
[TFHEの解説](https://nindanaoto.github.io/)に基づいたC++20によるTFHEライブラリ。

## det-WFA

### Build

On Ubuntu 20.04.2 LTS:

```sh
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_COMPILER="g++-10" ..
$ make
```

### Run

In the `build/` directory:

```sh
$ time ./main ../test/01.spec ../test/01.in
Parameter:
	Input size:	264
	State size:	19561
	Weight Num Scale:	1
	Weight size:	19561
	Concurrency:	12
=====
[263] #CMUX : 19537
[262] #CMUX : 19543
[261] #CMUX : 19537
[260] #CMUX : 19543
[259] #CMUX : 19537

...

[7] #CMUX : 81
[6] #CMUX : 62
[5] #CMUX : 39
[4] #CMUX : 21
[3] #CMUX : 10
[2] #CMUX : 5
[1] #CMUX : 2
[0] #CMUX : 1
Total #CMUX : 4483779
Result (bitstr):

00000000 00000000 00000000 00000000
00000001 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000001 00000000 00000000
00000000 00000000 00000000 00000001
00000000 00000000 00000000 00000001
00000000 00001000 00000001 00000000
00000000 00000000 00000000 00000001
00000001 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000
00000000 00000000 00000000 00000000

./main ../test/01.spec ../test/01.in  160.97s user 0.59s system 1054% cpu 15.317 total
```

The result means:
```
quick fox jumps over the lazy dog
000010000000010000010001001000011
```

## 章立てに対応するコミット

- [0.Introduction](https://github.com/ushitora-anqou/aqtfhe3/commit/37c416431c277d83341bc00a6ef8c62a588c76b2)
- [1.TLWE](https://github.com/ushitora-anqou/aqtfhe3/commit/41b6267cf46b9b616cc2bd718cddf96c7e5fe705)
- [2.TRLWE and Sample ExtractIndex](https://github.com/ushitora-anqou/aqtfhe3/commit/26684dbd50594370a257f192965b3111d9c9eb4c)
- [3.TRGSW and CMUX](https://github.com/ushitora-anqou/aqtfhe3/commit/8fcc5d9e0ec5fecc15abc90da5775727b2cd3f27)
- [4.Blind Rotate](https://github.com/ushitora-anqou/aqtfhe3/commit/a640edd78664286713f21c218e30fc74873532d2)
- [5.Identity Key Switch](https://github.com/ushitora-anqou/aqtfhe3/commit/c3e499b85a6c06f6f1604d3fdad40d9d39d143e1)
- [6.HomNAND](https://github.com/ushitora-anqou/aqtfhe3/commit/dc7cc98ca0b4d2be36a9d4eab45a25a669f50b2f)
- [7.Polynomial Multiplication By FFT](https://github.com/ushitora-anqou/aqtfhe3/commit/dacc5a7b0a4b93e52370d3519c305ad38a968735)
    - FFTには[GamaらのTFHEライブラリ](https://tfhe.github.io/tfhe/)に付属するspqliosを[Xbyak](https://github.com/herumi/xbyak)を用いてC++上に移植したものを使用しています。JITであることを活かした最適化は[今後の課題](https://github.com/ushitora-anqou/aqtfhe3/pull/1)です。

これらのコミットは、開発を終えた後に新しくgit treeを作成し整理したものです。
関係のないコード片が混ざっているコミットがあった場合は、
IssueやPull Requestなどでお知らせください。

## ビルド

```sh
$ make # ビルド
```

Ubuntu 20.04 LTSにおいて動作確認しています。
C++20に対応したコンパイラ（GCC 10やClang 10など）が必要です。
必要に応じて`Makefile`内で起動されるC++コンパイラを変更してください。

## 実行

```sh
$ make test # テスト及びベンチマーク
```

デフォルトでは[128-bit securityのパラメタ[CGGI19]](https://tfhe.github.io/tfhe/security_and_params.html)を使用します。
Intel i7-8700（最大周波数4.6GHz）で実行すると、例えば次のように出力されます。

```
Generating secret key (0) (SEED: 659303009)
==============================
0 (SEED: 1507937980)
test_tlwe:
	986510248 == 986510248
	0 == 0
test_trlwe:
	[0] 0 == 0
	[0] 0 == 0
test_cmux:
	[0] 0 ? 0 : 1 == 1
test_blind_rotate:
	0 == 0
test_identity_key_switch:
	1 == 1
test_hom_nand:
	nand(1, 0) = 1
==============================

Test passed. (1464 ms)
Benchmark started...
done.
Benchmark result: 12.642 ms/gate
```

## 書いた人

[Ushitora Anqou](https://anqou.net/)

## ライセンス

LICENSEファイルを参照してください。
