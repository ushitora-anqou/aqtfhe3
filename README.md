# aqTFHE3

[@nindanaoto](https://github.com/nindanaoto)による
[TFHEの解説](https://nindanaoto.github.io/)に基づいたC++20によるTFHEライブラリ。

## 章立てに対応するコミット

- [0.Introduction](https://github.com/ushitora-anqou/aqtfhe3/commit/37c416431c277d83341bc00a6ef8c62a588c76b2)
- [1.TLWE](https://github.com/ushitora-anqou/aqtfhe3/commit/41b6267cf46b9b616cc2bd718cddf96c7e5fe705)
- [2.TRLWE and Sample ExtractIndex](https://github.com/ushitora-anqou/aqtfhe3/commit/26684dbd50594370a257f192965b3111d9c9eb4c)
- [3.TRGSW and CMUX](https://github.com/ushitora-anqou/aqtfhe3/commit/36797f4744b168cf1a0ab09962d2b5151a49d169)
- [4.Blind Rotate](https://github.com/ushitora-anqou/aqtfhe3/commit/f3015dd673aee47bcefaa0b9ce2b8fb5a8385e9e)
- [5.Identity Key Switch](https://github.com/ushitora-anqou/aqtfhe3/commit/952c2c507a6f786493b6140786504ad6f93c1b84)
- [6.HomNAND](https://github.com/ushitora-anqou/aqtfhe3/commit/cde50713a72ffaf4d6b6797137915407ae007175)
- [7.Polynomial Multiplication By FFT](https://github.com/ushitora-anqou/aqtfhe3/commit/5a454eda735cc9e678e53f1be6ad6efdb5f48382)
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
