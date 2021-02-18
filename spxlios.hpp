/*
    spxlios.hpp is strongly inspired by spqlios processor in library
    TFHE (https://tfhe.github.io/tfhe/). See LICENSE file for its license.
*/
#pragma once
#ifndef SPXLIOS_HPP
#define SPXLIOS_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>

#include "xbyak/xbyak/xbyak.h"

namespace spxlios {

namespace detail {
struct alignas(32) FFT_PRECOMP {
    uint64_t n;
    double *trig_tables;
};

enum UnrollOption {
    UNROLL_FFT_1 = 1 << 0,
    UNROLL_FFT_2 = 1 << 1,
    UNROLL_FFT_3 = 1 << 2,
    UNROLL_FFT_4 = 1 << 3,
    UNROLL_FFT_5 = 1 << 4,

    UNROLL_IFFT_1 = 1 << 0,
    UNROLL_IFFT_2 = 1 << 1,
    UNROLL_IFFT_3 = 1 << 2,
    UNROLL_IFFT_4 = 1 << 3,
    UNROLL_IFFT_5 = 1 << 4,
};

struct fft_code : Xbyak::CodeGenerator {
    fft_code(size_t N, int unroll = UNROLL_FFT_1 | UNROLL_FFT_5)
        : Xbyak::CodeGenerator(4096, Xbyak::AutoGrow)
    {
        /*
           index  0    1    2    3    4      5
           reg    rdi, rsi, rdx, rcx, r8,    r9
           var    dst  sit  send bla, table, inout
        */
        using namespace Xbyak;
        const size_t n = 2 * N, ns4 = n / 4;

        vbroadcastsd(ymm2, ptr[rcx]);
        if (unroll & UNROLL_FFT_1) {
            // L("@@");
            for (size_t i = 0; i < N; i += 4) {
                vmovupd(ymm0, ptr[rsi]);
                vmulpd(ymm0, ymm2, ymm0);
                vmovapd(ptr[rdi], ymm0);
                add(rsi, 32);
                add(rdi, 32);
                // cmp(rsi, rdx);
                // jb("@b");
            }
        }
        else {
            L("@@");
            vmovupd(ymm0, ptr[rsi]);
            vmulpd(ymm0, ymm2, ymm0);
            vmovapd(ptr[rdi], ymm0);
            add(rsi, 32);
            add(rdi, 32);
            cmp(rsi, rdx);
            jb("@b");
        }

        mov(rdi, r8);
        mov(rsi, r9);

        push(r10);
        push(r11);
        push(r12);
        push(r13);
        push(r14);
        push(rbx);

        mov(rax, rdi);
        mov(rdi, rsi);

        mov(rdx, qword[rax]);
        mov(r8, qword[rax + 0x8]);

        mov(r9, rdx);
        shr(r9, 0x2);
        lea(rsi, ptr[rdi + r9 * 8]);

        Label size4negation0, size4negation1, size4negation2, size4negation3;
        vmovapd(ymm15, ptr[rip + size4negation0]);
        vmovapd(ymm14, ptr[rip + size4negation1]);
        vmovapd(ymm13, ptr[rip + size4negation2]);
        vmovapd(ymm12, ptr[rip + size4negation3]);

        if (unroll & UNROLL_FFT_2) {
            // mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            // L("fftsize2loop");
            for (size_t block = 0; block < ns4; block += 4) {
                vmovapd(ymm0, ptr[r10]);
                vmovapd(ymm1, ptr[r11]);
                vshufpd(ymm2, ymm0, ymm0, 0x0);
                vshufpd(ymm3, ymm0, ymm0, 0xf);
                vshufpd(ymm4, ymm1, ymm1, 0x0);
                vshufpd(ymm5, ymm1, ymm1, 0xf);
                vfmadd231pd(ymm2, ymm12, ymm3);
                vfmadd231pd(ymm4, ymm12, ymm5);
                vmovapd(ptr[r10], ymm2);
                vmovapd(ptr[r11], ymm4);
                lea(r10, ptr[r10 + 0x20]);
                lea(r11, ptr[r11 + 0x20]);
                // add(rax, 0x4);
                // cmp(rax, r9);
                // jb("fftsize2loop");
            }
        }
        else {
            mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            L("fftsize2loop");
            vmovapd(ymm0, ptr[r10]);
            vmovapd(ymm1, ptr[r11]);
            vshufpd(ymm2, ymm0, ymm0, 0x0);
            vshufpd(ymm3, ymm0, ymm0, 0xf);
            vshufpd(ymm4, ymm1, ymm1, 0x0);
            vshufpd(ymm5, ymm1, ymm1, 0xf);
            vfmadd231pd(ymm2, ymm12, ymm3);
            vfmadd231pd(ymm4, ymm12, ymm5);
            vmovapd(ptr[r10], ymm2);
            vmovapd(ptr[r11], ymm4);
            lea(r10, ptr[r10 + 0x20]);
            lea(r11, ptr[r11 + 0x20]);
            add(rax, 0x4);
            cmp(rax, r9);
            jb("fftsize2loop");
        }

        if (unroll & UNROLL_FFT_3) {
            // mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            // L("fftsize4loop");
            for (size_t block = 0; block < ns4; block += 4) {
                vmovapd(ymm0, ptr[r10]);
                vmovapd(ymm1, ptr[r11]);
                vperm2f128(ymm4, ymm0, ymm0, 0x20);
                vperm2f128(ymm5, ymm1, ymm1, 0x20);
                vperm2f128(ymm6, ymm0, ymm0, 0x31);
                vperm2f128(ymm7, ymm1, ymm1, 0x31);
                vshufpd(ymm8, ymm6, ymm7, 0xa);
                vshufpd(ymm9, ymm7, ymm6, 0xa);
                vfmadd231pd(ymm4, ymm13, ymm8);
                vfmadd231pd(ymm5, ymm14, ymm9);
                vmovapd(ptr[r10], ymm4);
                vmovapd(ptr[r11], ymm5);
                lea(r10, ptr[r10 + 0x20]);
                lea(r11, ptr[r11 + 0x20]);
                // add(rax, 0x4);
                // cmp(rax, r9);
                // jb("fftsize4loop");
            }
        }
        else {
            mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            L("fftsize4loop");
            vmovapd(ymm0, ptr[r10]);
            vmovapd(ymm1, ptr[r11]);
            vperm2f128(ymm4, ymm0, ymm0, 0x20);
            vperm2f128(ymm5, ymm1, ymm1, 0x20);
            vperm2f128(ymm6, ymm0, ymm0, 0x31);
            vperm2f128(ymm7, ymm1, ymm1, 0x31);
            vshufpd(ymm8, ymm6, ymm7, 0xa);
            vshufpd(ymm9, ymm7, ymm6, 0xa);
            vfmadd231pd(ymm4, ymm13, ymm8);
            vfmadd231pd(ymm5, ymm14, ymm9);
            vmovapd(ptr[r10], ymm4);
            vmovapd(ptr[r11], ymm5);
            lea(r10, ptr[r10 + 0x20]);
            lea(r11, ptr[r11 + 0x20]);
            add(rax, 0x4);
            cmp(rax, r9);
            jb("fftsize4loop");
        }

        if (unroll & UNROLL_FFT_4) {
            mov(rdx, r8);
            mov(rax, 0x4);
            // L("ffthalfnnloop");
            for (size_t halfnn = 4; halfnn < ns4; halfnn *= 2) {
                inLocalLabel();
                // size_t nn = 2 * halfnn;
                mov(rbx, 0x0);
                L(".fftblockloop");
                // for (size_t block = 0; block < ns4; block += nn) {
                lea(r10, ptr[rdi + rbx * 8]);
                lea(r11, ptr[rsi + rbx * 8]);
                lea(r12, ptr[r10 + rax * 8]);
                lea(r13, ptr[r11 + rax * 8]);
                mov(r14, rdx);
                mov(rcx, 0x0);
                L(".fftoffloop");
                // for (size_t off = 0; off < halfnn; off += 4) {
                vmovapd(ymm0, ptr[r10]);
                vmovapd(ymm1, ptr[r11]);
                vmovapd(ymm2, ptr[r12]);
                vmovapd(ymm3, ptr[r13 + 0x0]);
                vmovapd(ymm4, ptr[r14]);
                vmovapd(ymm5, ptr[r14 + 0x20]);
                vmulpd(ymm6, ymm4, ymm2);
                vmulpd(ymm7, ymm5, ymm2);
                vfnmadd231pd(ymm6, ymm5, ymm3);
                vfmadd231pd(ymm7, ymm4, ymm3);
                vsubpd(ymm2, ymm0, ymm6);
                vsubpd(ymm3, ymm1, ymm7);
                vaddpd(ymm0, ymm0, ymm6);
                vaddpd(ymm1, ymm1, ymm7);
                vmovapd(ptr[r10], ymm0);
                vmovapd(ptr[r11], ymm1);
                vmovapd(ptr[r12], ymm2);
                vmovapd(ptr[r13 + 0x0], ymm3);
                lea(r10, ptr[r10 + 0x20]);
                lea(r11, ptr[r11 + 0x20]);
                lea(r12, ptr[r12 + 0x20]);
                lea(r13, ptr[r13 + 0x20]);
                lea(r14, ptr[r14 + 0x40]);
                add(rcx, 0x4);
                cmp(rcx, rax);
                jb(".fftoffloop");
                //}
                lea(rbx, ptr[rbx + rax * 2]);
                cmp(rbx, r9);
                jb(".fftblockloop");
                //}
                shl(rax, 1);
                lea(rdx, ptr[rdx + rax * 8]);
                // cmp(rax, r9);
                // jb("ffthalfnnloop");
                outLocalLabel();
            }
        }
        else {
            mov(rdx, r8);
            mov(rax, 0x4);
            L("ffthalfnnloop");
            mov(rbx, 0x0);
            L("fftblockloop");
            lea(r10, ptr[rdi + rbx * 8]);
            lea(r11, ptr[rsi + rbx * 8]);
            lea(r12, ptr[r10 + rax * 8]);
            lea(r13, ptr[r11 + rax * 8]);
            mov(r14, rdx);
            mov(rcx, 0x0);
            L("fftoffloop");
            vmovapd(ymm0, ptr[r10]);
            vmovapd(ymm1, ptr[r11]);
            vmovapd(ymm2, ptr[r12]);
            vmovapd(ymm3, ptr[r13 + 0x0]);
            vmovapd(ymm4, ptr[r14]);
            vmovapd(ymm5, ptr[r14 + 0x20]);
            vmulpd(ymm6, ymm4, ymm2);
            vmulpd(ymm7, ymm5, ymm2);
            vfnmadd231pd(ymm6, ymm5, ymm3);
            vfmadd231pd(ymm7, ymm4, ymm3);
            vsubpd(ymm2, ymm0, ymm6);
            vsubpd(ymm3, ymm1, ymm7);
            vaddpd(ymm0, ymm0, ymm6);
            vaddpd(ymm1, ymm1, ymm7);
            vmovapd(ptr[r10], ymm0);
            vmovapd(ptr[r11], ymm1);
            vmovapd(ptr[r12], ymm2);
            vmovapd(ptr[r13 + 0x0], ymm3);
            lea(r10, ptr[r10 + 0x20]);
            lea(r11, ptr[r11 + 0x20]);
            lea(r12, ptr[r12 + 0x20]);
            lea(r13, ptr[r13 + 0x20]);
            lea(r14, ptr[r14 + 0x40]);
            add(rcx, 0x4);
            cmp(rcx, rax);
            jb("fftoffloop");
            lea(rbx, ptr[rbx + rax * 2]);
            cmp(rbx, r9);
            jb("fftblockloop");
            shl(rax, 1);
            lea(rdx, ptr[rdx + rax * 8]);
            cmp(rax, r9);
            jb("ffthalfnnloop");
        }

        if (unroll & UNROLL_FFT_5) {
            // mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            // L("fftfinalloop");
            for (size_t j = 0; j < ns4; j += 4) {
                vmovapd(ymm0, ptr[r10]);
                vmovapd(ymm1, ptr[r11]);
                vmovapd(ymm2, ptr[rdx]);
                vmovapd(ymm3, ptr[rdx + 0x20]);
                vmulpd(ymm4, ymm2, ymm0);
                vmulpd(ymm5, ymm3, ymm0);
                vmulpd(ymm6, ymm2, ymm1);
                vmulpd(ymm7, ymm3, ymm1);
                vsubpd(ymm0, ymm4, ymm7);
                vaddpd(ymm1, ymm5, ymm6);
                vmovapd(ptr[r10], ymm0);
                vmovapd(ptr[r11], ymm1);
                lea(r10, ptr[r10 + 0x20]);
                lea(r11, ptr[r11 + 0x20]);
                lea(rdx, ptr[rdx + 0x40]);
                // add(rax, 0x4);
                // cmp(rax, r9);
                // jb("fftfinalloop");
            }
        }
        else {
            mov(rax, 0x0);
            mov(r10, rdi);
            mov(r11, rsi);
            L("fftfinalloop");
            vmovapd(ymm0, ptr[r10]);
            vmovapd(ymm1, ptr[r11]);
            vmovapd(ymm2, ptr[rdx]);
            vmovapd(ymm3, ptr[rdx + 0x20]);
            vmulpd(ymm4, ymm2, ymm0);
            vmulpd(ymm5, ymm3, ymm0);
            vmulpd(ymm6, ymm2, ymm1);
            vmulpd(ymm7, ymm3, ymm1);
            vsubpd(ymm0, ymm4, ymm7);
            vaddpd(ymm1, ymm5, ymm6);
            vmovapd(ptr[r10], ymm0);
            vmovapd(ptr[r11], ymm1);
            lea(r10, ptr[r10 + 0x20]);
            lea(r11, ptr[r11 + 0x20]);
            lea(rdx, ptr[rdx + 0x40]);
            add(rax, 0x4);
            cmp(rax, r9);
            jb("fftfinalloop");
        }

        // L("fftend");
        vzeroall();
        pop(rbx);
        pop(r14);
        pop(r13);
        pop(r12);
        pop(r11);
        pop(r10);
        ret();

        align(32);

        // size4negation0: .double +1.0, +1.0, +1.0, -1.0 /* ymm15 */
        L(size4negation0);
        const double tbl0[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  //-1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl0[i]);

        // size4negation1: .double +1.0, -1.0, -1.0, +1.0 /* ymm14 */
        L(size4negation1);
        const double tbl1[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl1[i]);

        // size4negation2: .double +1.0, +1.0, -1.0, -1.0 /* ymm13 */
        L(size4negation2);
        const double tbl2[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl2[i]);

        // size4negation3: .double +1.0, -1.0, +1.0, -1.0 /* ymm12 */
        L(size4negation3);
        const double tbl3[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl3[i]);
    }
};

struct ifft_code : Xbyak::CodeGenerator {
    ifft_code(size_t N, int unroll = UNROLL_IFFT_1 | UNROLL_IFFT_5)
        : Xbyak::CodeGenerator(4096, Xbyak::AutoGrow)
    {
        /*
           index  0    1    2     3      4      5   | 6,   7,   8
           reg    rdi, rsi, rdx,  rcx,   r8,    r9  |
           var    dst  ait, aend, table, inout, _,  | dst, sit, send
        */
        using namespace Xbyak;
        const size_t n = 2 * N, ns4 = n / 4;

        if (unroll & UNROLL_IFFT_1) {
            // L("@@");
            for (size_t i = 0; i < N; i += 4) {
                vmovupd(xmm0, ptr[rsi]);
                vcvtdq2pd(ymm1, xmm0);
                vmovapd(ptr[rdi], ymm1);
                add(rsi, 16);
                add(rdi, 32);
                // cmp(rsi, rdx);
                // jb("@b");
            }
        }
        else {
            L("@@");
            vmovupd(xmm0, ptr[rsi]);
            vcvtdq2pd(ymm1, xmm0);
            vmovapd(ptr[rdi], ymm1);
            add(rsi, 16);
            add(rdi, 32);
            cmp(rsi, rdx);
            jb("@b");
        }

        mov(rdi, rcx);
        mov(rsi, r8);

        L("ifft");
        push(r10);
        push(r11);
        push(r12);
        push(r13);
        push(r14);
        push(rbx);

        mov(rax, rdi);
        mov(rdi, rsi);

        mov(rdx, qword[rax]);
        mov(r8, qword[rax + 0x8]);

        mov(r10, rdx);
        shl(r10, 1);
        add(rsi, r10);

        if (unroll & UNROLL_IFFT_2) {
            shr(r10, 0x3);
            mov(rcx, 0x0);
            mov(r11, r8);
            // L("firstloop");
            for (size_t j = 0; j < ns4; j += 4) {
                vmovapd(ymm0, ptr[rdi + rcx * 8]);
                vmovapd(ymm1, ptr[rsi + rcx * 8]);
                vmovapd(ymm2, ptr[r11]);
                vmovapd(ymm3, ptr[r11 + 0x20]);
                vmulpd(ymm4, ymm2, ymm0);
                vmulpd(ymm5, ymm3, ymm0);
                vfnmadd231pd(ymm4, ymm3, ymm1);
                vfmadd231pd(ymm5, ymm2, ymm1);
                vmovapd(ptr[rdi + rcx * 8], ymm4);
                vmovapd(ptr[rsi + rcx * 8], ymm5);
                lea(r11, ptr[r11 + 0x40]);
                add(rcx, 0x4);
                // cmp(rcx, r10);
                // jb("firstloop");
            }
        }
        else {
            shr(r10, 0x3);
            mov(rcx, 0x0);
            mov(r11, r8);
            L("firstloop");
            vmovapd(ymm0, ptr[rdi + rcx * 8]);
            vmovapd(ymm1, ptr[rsi + rcx * 8]);
            vmovapd(ymm2, ptr[r11]);
            vmovapd(ymm3, ptr[r11 + 0x20]);
            vmulpd(ymm4, ymm2, ymm0);
            vmulpd(ymm5, ymm3, ymm0);
            vfnmadd231pd(ymm4, ymm3, ymm1);
            vfmadd231pd(ymm5, ymm2, ymm1);
            vmovapd(ptr[rdi + rcx * 8], ymm4);
            vmovapd(ptr[rsi + rcx * 8], ymm5);
            lea(r11, ptr[r11 + 0x40]);
            add(rcx, 0x4);
            cmp(rcx, r10);
            jb("firstloop");
        }

        if (unroll & UNROLL_IFFT_3) {
            mov(r12, r10);
            // L("nnloop");
            for (size_t nn = ns4; nn >= 8; nn /= 2) {
                inLocalLabel();
                // size_t halfnn = nn / 2;
                mov(r13, r12);
                shr(r13, 1);
                lea(r8, ptr[r8 + r12 * 8]);
                lea(r8, ptr[r8 + r12 * 8]);
                mov(r11, 0x0);
                L(".blockloop");
                // for (size_t block = 0; block < ns4; block += nn) {
                lea(rax, ptr[rdi + r11 * 8]);
                lea(rbx, ptr[rsi + r11 * 8]);
                lea(rcx, ptr[rax + r13 * 8]);
                lea(rdx, ptr[rbx + r13 * 8]);
                mov(r9, 0x0);
                mov(r14, r8);
                L(".offloop");
                // for (size_t off = 0; off < halfnn; off += 4) {
                vmovapd(ymm0, ptr[rax + r9 * 8]);
                vmovapd(ymm1, ptr[rbx + r9 * 8]);
                vmovapd(ymm2, ptr[rcx + r9 * 8]);
                vmovapd(ymm3, ptr[rdx + r9 * 8]);
                vaddpd(ymm4, ymm2, ymm0);
                vaddpd(ymm5, ymm3, ymm1);
                vsubpd(ymm6, ymm0, ymm2);
                vsubpd(ymm7, ymm1, ymm3);
                vmovapd(ptr[rax + r9 * 8], ymm4);
                vmovapd(ptr[rbx + r9 * 8], ymm5);
                vmovapd(ymm8, ptr[r14]);
                vmovapd(ymm9, ptr[r14 + 0x20]);
                vmulpd(ymm4, ymm8, ymm6);
                vfnmadd231pd(ymm4, ymm9, ymm7);
                vmulpd(ymm5, ymm9, ymm6);
                vfmadd231pd(ymm5, ymm8, ymm7);
                vmovapd(ptr[rcx + r9 * 8], ymm4);
                vmovapd(ptr[rdx + r9 * 8], ymm5);
                lea(r14, ptr[r14 + 0x40]);
                add(r9, 0x4);
                cmp(r9, r13);
                jb(".offloop");
                //}
                add(r11, r12);
                cmp(r11, r10);
                jb(".blockloop");
                //}
                mov(r12, r13);
                // cmp(r12, 0x8);
                outLocalLabel();
            }
        }
        else {
            mov(r12, r10);
            L("nnloop");
            mov(r13, r12);
            shr(r13, 1);
            lea(r8, ptr[r8 + r12 * 8]);
            lea(r8, ptr[r8 + r12 * 8]);
            mov(r11, 0x0);
            L("blockloop");
            lea(rax, ptr[rdi + r11 * 8]);
            lea(rbx, ptr[rsi + r11 * 8]);
            lea(rcx, ptr[rax + r13 * 8]);
            lea(rdx, ptr[rbx + r13 * 8]);
            mov(r9, 0x0);
            mov(r14, r8);
            L("offloop");
            vmovapd(ymm0, ptr[rax + r9 * 8]);
            vmovapd(ymm1, ptr[rbx + r9 * 8]);
            vmovapd(ymm2, ptr[rcx + r9 * 8]);
            vmovapd(ymm3, ptr[rdx + r9 * 8]);
            vaddpd(ymm4, ymm2, ymm0);
            vaddpd(ymm5, ymm3, ymm1);
            vsubpd(ymm6, ymm0, ymm2);
            vsubpd(ymm7, ymm1, ymm3);
            vmovapd(ptr[rax + r9 * 8], ymm4);
            vmovapd(ptr[rbx + r9 * 8], ymm5);
            vmovapd(ymm8, ptr[r14]);
            vmovapd(ymm9, ptr[r14 + 0x20]);
            vmulpd(ymm4, ymm8, ymm6);
            vfnmadd231pd(ymm4, ymm9, ymm7);
            vmulpd(ymm5, ymm9, ymm6);
            vfmadd231pd(ymm5, ymm8, ymm7);
            vmovapd(ptr[rcx + r9 * 8], ymm4);
            vmovapd(ptr[rdx + r9 * 8], ymm5);
            lea(r14, ptr[r14 + 0x40]);
            add(r9, 0x4);
            cmp(r9, r13);
            jb("offloop");
            add(r11, r12);
            cmp(r11, r10);
            jb("blockloop");
            mov(r12, r13);
            cmp(r12, 0x8);
            jae("nnloop");
        }

        Label size4negation0, size4negation1, size4negation2, size4negation3;
        vmovapd(ymm15, ptr[rip + size4negation0]);
        vmovapd(ymm14, ptr[rip + size4negation1]);
        vmovapd(ymm13, ptr[rip + size4negation2]);
        vmovapd(ymm12, ptr[rip + size4negation3]);

        if (unroll & UNROLL_IFFT_3) {
            // mov(rax, 0x0);
            mov(r11, rdi);
            mov(r12, rsi);
            // L("size4loop");
            for (size_t block = 0; block < ns4; block += 4) {
                vmovapd(ymm0, ptr[r11]);
                vmovapd(ymm1, ptr[r12]);
                vshufpd(ymm2, ymm0, ymm1, 0xa);
                vshufpd(ymm3, ymm1, ymm0, 0xa);
                vperm2f128(ymm4, ymm0, ymm2, 0x20);
                vperm2f128(ymm5, ymm0, ymm2, 0x31);
                vperm2f128(ymm6, ymm1, ymm3, 0x20);
                vperm2f128(ymm7, ymm1, ymm3, 0x31);
                vmulpd(ymm4, ymm15, ymm4);
                vfmadd231pd(ymm4, ymm14, ymm5);
                vfmadd231pd(ymm6, ymm13, ymm7);
                vmovapd(ptr[r11], ymm4);
                vmovapd(ptr[r12], ymm6);
                lea(r11, ptr[r11 + 0x20]);
                lea(r12, ptr[r12 + 0x20]);
                // add(rax, 0x4);
                // cmp(rax, r10);
                // jb("size4loop");
            }
        }
        else {
            mov(rax, 0x0);
            mov(r11, rdi);
            mov(r12, rsi);
            L("size4loop");
            vmovapd(ymm0, ptr[r11]);
            vmovapd(ymm1, ptr[r12]);
            vshufpd(ymm2, ymm0, ymm1, 0xa);
            vshufpd(ymm3, ymm1, ymm0, 0xa);
            vperm2f128(ymm4, ymm0, ymm2, 0x20);
            vperm2f128(ymm5, ymm0, ymm2, 0x31);
            vperm2f128(ymm6, ymm1, ymm3, 0x20);
            vperm2f128(ymm7, ymm1, ymm3, 0x31);
            vmulpd(ymm4, ymm15, ymm4);
            vfmadd231pd(ymm4, ymm14, ymm5);
            vfmadd231pd(ymm6, ymm13, ymm7);
            vmovapd(ptr[r11], ymm4);
            vmovapd(ptr[r12], ymm6);
            lea(r11, ptr[r11 + 0x20]);
            lea(r12, ptr[r12 + 0x20]);
            add(rax, 0x4);
            cmp(rax, r10);
            jb("size4loop");
        }

        if (unroll & UNROLL_IFFT_4) {
            // mov(rax, 0x0);
            mov(r11, rdi);
            mov(r12, rsi);
            // L("size2loop");
            for (size_t block = 0; block < ns4; block += 4) {
                vmovapd(ymm0, ptr[r11]);
                vmovapd(ymm1, ptr[r12]);
                vshufpd(ymm2, ymm0, ymm0, 0x0);
                vshufpd(ymm3, ymm0, ymm0, 0xf);
                vshufpd(ymm4, ymm1, ymm1, 0x0);
                vshufpd(ymm5, ymm1, ymm1, 0xf);
                vfmadd231pd(ymm2, ymm12, ymm3);
                vfmadd231pd(ymm4, ymm12, ymm5);
                vmovapd(ptr[r11], ymm2);
                vmovapd(ptr[r12], ymm4);
                lea(r11, ptr[r11 + 0x20]);
                lea(r12, ptr[r12 + 0x20]);
                // add(rax, 0x4);
                // cmp(rax, r10);
                // jb("size2loop");
            }
        }
        else {
            mov(rax, 0x0);
            mov(r11, rdi);
            mov(r12, rsi);
            L("size2loop");
            vmovapd(ymm0, ptr[r11]);
            vmovapd(ymm1, ptr[r12]);
            vshufpd(ymm2, ymm0, ymm0, 0x0);
            vshufpd(ymm3, ymm0, ymm0, 0xf);
            vshufpd(ymm4, ymm1, ymm1, 0x0);
            vshufpd(ymm5, ymm1, ymm1, 0xf);
            vfmadd231pd(ymm2, ymm12, ymm3);
            vfmadd231pd(ymm4, ymm12, ymm5);
            vmovapd(ptr[r11], ymm2);
            vmovapd(ptr[r12], ymm4);
            lea(r11, ptr[r11 + 0x20]);
            lea(r12, ptr[r12 + 0x20]);
            add(rax, 0x4);
            cmp(rax, r10);
            jb("size2loop");
        }

        L("end");
        vzeroall();
        pop(rbx);
        pop(r14);
        pop(r13);
        pop(r12);
        pop(r11);
        pop(r10);

        mov(rdi, ptr[rsp + 8]);   // dst
        mov(rsi, ptr[rsp + 16]);  // sit
        mov(rdx, ptr[rsp + 24]);  // send

        if (unroll & UNROLL_IFFT_5) {
            // L("@@");
            for (size_t i = 0; i < N; i += 4) {
                vmovapd(ymm0, ptr[rsi]);
                vmovupd(ptr[rdi], ymm0);
                add(rsi, 32);
                add(rdi, 32);
                // cmp(rsi, rdx);
                // jb("@b");
            }
        }
        else {
            L("@@");
            vmovapd(ymm0, ptr[rsi]);
            vmovupd(ptr[rdi], ymm0);
            add(rsi, 32);
            add(rdi, 32);
            cmp(rsi, rdx);
            jb("@b");
        }
        vzeroall();

        ret();

        align(32);

        // size4negation0: .double +1.0, +1.0, +1.0, -1.0
        L(size4negation0);
        const double tbl0[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  //+1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  //-1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl0[i]);

        // size4negation1: .double +1.0, +1.0, -1.0, +1.0
        L(size4negation1);
        const double tbl1[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl1[i]);

        // size4negation2: .double +1.0, +1.0, -1.0, -1.0
        L(size4negation2);
        const double tbl2[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl2[i]);

        // size4negation3: .double +1.0, -1.0, +1.0, -1.0
        L(size4negation3);
        const double tbl3[] = {
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  // +1.0
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf,  // -1.0
        };
        for (size_t i = 0; i < 32; i++)
            db(tbl3[i]);
    }
};

double accurate_cos(int32_t i, int32_t n)  // cos(2pi*i/n)
{
    constexpr double pi = std::numbers::pi;
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4)
        return std::cos(2. * pi * (n - i) / static_cast<double>(n));
    if (i >= 2 * n / 4)
        return -std::cos(2. * pi * (i - n / 2) / static_cast<double>(n));
    if (i >= 1 * n / 4)
        return -std::cos(2. * pi * (n / 2 - i) / static_cast<double>(n));
    return std::cos(2. * pi * (i) / static_cast<double>(n));
}

double accurate_sin(int32_t i, int32_t n)  // sin(2pi*i/n)
{
    constexpr double pi = std::numbers::pi;
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4)
        return -std::sin(2. * pi * (n - i) / static_cast<double>(n));
    if (i >= 2 * n / 4)
        return -std::sin(2. * pi * (i - n / 2) / static_cast<double>(n));
    if (i >= 1 * n / 4)
        return std::sin(2. * pi * (n / 2 - i) / static_cast<double>(n));
    return std::sin(2. * pi * (i) / static_cast<double>(n));
}

template <size_t N>
double *direct_trig_tables()
{
    constexpr int32_t n = 2 * N, ns4 = n / 4;
    double *head = new (std::align_val_t(32)) double[n];
    double *ptr = head;

    // subsequent iterations
    for (int32_t halfnn = 4; halfnn < ns4; halfnn *= 2) {
        int32_t nn = 2 * halfnn;
        int32_t j = n / nn;
        // cerr << "- b: " << halfnn  << "(offset: " <<
        // (ptr-reps->trig_tables) << ", mult: " << j << ")" << endl;
        for (int32_t i = 0; i < halfnn; i += 4) {
            // cerr << "--- i: " << i << endl;
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = detail::accurate_cos(-j * (i + k), n);
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = detail::accurate_sin(-j * (i + k), n);
        }
    }
    // last iteration
    for (int32_t i = 0; i < ns4; i += 4) {
        for (int32_t k = 0; k < 4; k++)
            *(ptr++) = detail::accurate_cos(-(i + k), n);
        for (int32_t k = 0; k < 4; k++)
            *(ptr++) = detail::accurate_sin(-(i + k), n);
    }

    return head;
}

template <size_t N>
double *reverse_trig_tables()
{
    constexpr int32_t n = 2 * N, ns4 = n / 4;
    double *head = new (std::align_val_t(32)) double[n];
    double *ptr = head;

    // first iteration
    for (int32_t j = 0; j < ns4; j += 4) {
        for (int32_t k = 0; k < 4; k++)
            *(ptr++) = detail::accurate_cos(j + k, n);
        for (int32_t k = 0; k < 4; k++)
            *(ptr++) = detail::accurate_sin(j + k, n);
    }
    // subsequent iterations
    for (int32_t nn = ns4; nn >= 8; nn /= 2) {
        int32_t halfnn = nn / 2;
        int32_t j = n / nn;
        // cerr << "- b: " << nn  << "(offset: " <<
        // (ptr-reps->trig_tables) << ", mult: " << j << ")" << endl;
        for (int32_t i = 0; i < halfnn; i += 4) {
            // cerr << "--- i: " << i << endl;
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = detail::accurate_cos(j * (i + k), n);
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = detail::accurate_sin(j * (i + k), n);
        }
    }

    return head;
}

}  // namespace detail

template <size_t N>
class processor {
    static_assert(N >= 16, "N must be >=16");
    static_assert((N & (N - 1)) == 0, "N must be a power of 2");

private:
    std::shared_ptr<double[]> direct_trig_tables_, reverse_trig_tables_;
    detail::FFT_PRECOMP tables_direct_, tables_reverse_;
    detail::fft_code fft_code_;
    detail::ifft_code ifft_code_;
    void (*fft_)(double *, const double *, const double *, const double *,
                 const detail::FFT_PRECOMP *, double *);
    void (*ifft_)(double *, const int32_t *, const int32_t *,
                  const detail::FFT_PRECOMP *, double *, int, double *,
                  double *, double *);

public:
    processor()
        : direct_trig_tables_(detail::direct_trig_tables<N>()),
          reverse_trig_tables_(detail::reverse_trig_tables<N>()),
          tables_direct_(),
          tables_reverse_(),
          fft_code_(N, detail::UNROLL_FFT_1),
          ifft_code_(N, detail::UNROLL_IFFT_1),
          fft_(nullptr),
          ifft_(nullptr)
    {
        tables_direct_.n = 2 * N;
        tables_direct_.trig_tables = direct_trig_tables_.get();
        tables_reverse_.n = 2 * N;
        tables_reverse_.trig_tables = reverse_trig_tables_.get();
        fft_code_.ready();
        fft_ = fft_code_.getCode<void (*)(
            double *, const double *, const double *, const double *,
            const detail::FFT_PRECOMP *, double *)>();
        ifft_code_.ready();
        ifft_ =
            ifft_code_
                .getCode<void (*)(double *, const int32_t *, const int32_t *,
                                  const detail::FFT_PRECOMP *, double *, int,
                                  double *, double *, double *)>();
    }

    void execute_direct_torus32(std::array<uint32_t, N> &out,
                                const std::array<double, N> &src)
    {
        uint32_t *res = out.data();
        const double *a = src.data();
        alignas(32) double real_inout_direct[N];

        constexpr double _2sN = 2.0 / N;
        double *dst = real_inout_direct;
        const double *sit = a;
        const double *send = a + N;
        const double *bla = &_2sN;

        fft_(dst, sit, send, bla, &tables_direct_, real_inout_direct);

        for (size_t i = 0; i < N; i++)
            res[i] = uint32_t(int64_t(real_inout_direct[i]));
    }

    void execute_reverse_torus32(std::array<double, N> &out,
                                 const std::array<uint32_t, N> &src)
    {
        double *res = out.data();
        const int32_t *a = reinterpret_cast<const int32_t *>(src.data());
        alignas(32) double real_inout_rev[N];

        double *dst = real_inout_rev;
        const int32_t *ait = a;
        const int32_t *aend = a + N;

        double *dst2 = res;
        double *sit = real_inout_rev;
        double *send = real_inout_rev + N;

        ifft_(dst, ait, aend, &tables_reverse_, real_inout_rev, 0, dst2, sit,
              send);
    }
};  // namespace spxlios

}  // namespace spxlios

#endif
