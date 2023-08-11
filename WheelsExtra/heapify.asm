

COMMENT @***
 - INVOCATION CODE -
 mov  r9,	qword ptr [pRoot]
 mov  r8,	[PrimeEngine::vSieve]  
 mov  edx,  dword ptr [PrimeEngine::idx_i]  
 mov  ecx,  dword ptr [PrimeEngine::idx_j]  
 call HeapifyASM 

 - INPUTS -
	EDX	=>	idx_i
	ECX	=>	idx_j
	R9	=>	pRoot
	R8	=>	vSieve

 - LOCALS -
	RSI	=>	vSieve
	RDI => 	pChildL or pChildR (after choosing the lowest)
	
	R8	=>	idxj_child
	R9	=>	idxi_child (after choosing the lowest)
	
	R10 =>  last_idx_j
	R11 =>  last_idx_i
	
	R12 =>  idx_j_current
	R13 =>  idx_i_current
	
	RBX =>	valnow

	R14	=>	idxj_child_R
	R15	=>	idxi_child_R

@***

.data

pRoot		DQ		0

.const
ROW_SIZE		EQU		1000000
LAST_ROW_IDX	EQU		ROW_SIZE-1

.code

HeapifyASM proc
	;push rax					
	push rbx		; x64 calling convention:
	;push rcx		; RAX, RCX, RDX, R8, R9, R10, R11 are considered volatile
	;push rdx		; RBX, RBP, RDI, RSI, RSP, R12, R13, R14, and R15 are considered nonvolatile
	push rdi
	push rsi
	;push r8
	;push r9
	;push r10
	;push r11
	push r12
	push r13
	push r14
	push r15
	push rbp

; initialize idx_current with {0, 0}
    xor		r12, r12
    xor		r13, r13

; initialize vSieve
	mov		rsi, r8
; store pRoot
	;mov pRoot, r9				
	mov rbp, r9				
; initialize valnow
	mov rbx, qword ptr [r9]		; CurrentValue is at offset 0
	
; initialize idx_child with {0, 1}
	xor 	r8, r8
	mov		r9, 1

; initialize last_idx_j
	mov		r10d, ecx
; initialize last_idx_i
	test    dl,1			; test if idx_i (edx) is EVEN
	lea     r11d,[edx-1]	; for EVEN case: load EDX-1
	cmove   r11d,edx		; but keep the original value if ODD
							; use cmove here because branching is unpredictable...

LBL_START_LOOP:

; vSieve[idxj_child_L] => rax (current row)
	mov     rax, qword ptr [rsi+r8*8]
; vSieve[idxj_child_L][idxi_child_L] => rdi (pChildL)
	mov     rdi, qword ptr [rax+r9*8]
; vSieve[idxj_child_L][idxi_child_L].uCurrentValue => rdx
	mov     rdx, qword ptr [rdi]
; now we have left child value in rdx and left child in rdi!

; if at the end of the row we must advance to the next one
	cmp		r9, LAST_ROW_IDX
	je		LBL_RIGHT_CHILD_ON_NEXT_ROW		; it can not be greater!!!
; jump if at the end of the row
; otherwise right child is the next in current row

; if here (more probable branch, although only marginal)
; set values for right child indexes
; idxj_child_R = idxj_child_L (same row)
	mov		r14, r8
; idxi_child_R = idxi_child_L + 1 (next child)
	lea		r15, [r9+1]
; vSieve[idxj_child_R] is already in rax (current row), from above (jR==jL)
; vSieve[idxj_child_R][idxi_child_R] => rax (pChildR)
	mov     rax, qword ptr [rax+r15*8]
; now we have right child in rax!
	
; 	jump over next session, to comparison
	jmp		LBL_COMPARE_LEFT_RIGHT	

LBL_RIGHT_CHILD_ON_NEXT_ROW:
; idxj_child_R = idxj_child_L+1
	lea		r14, [r8+1]
; idxi_child_R = 0
	xor		r15, r15
; vSieve[idxj_child_R] => rax
	mov     rax, qword ptr [rsi+r14*8]
; vSieve[idxj_child_R][idxi_child_R] => rax (pChildR)
	mov     rax, qword ptr [rax]	; r15 is 0 here
; now we have right child in rax!

LBL_COMPARE_LEFT_RIGHT:
; at this point we have L indexes in r8, 9 / R in r14, 15
; L in rdi / R in rax
; vSieve[idxj_child_R][idxi_child_R].uCurrentValue => rcx
	mov		rcx, qword ptr [rax]
; and values in rdx and rcx
	cmp rdx, rcx
	jle	LBL_COMPARE_WITH_valnow

;if here, R is smaller, so we go with R as current
	mov 	rdx, rcx 	; value
	mov 	rdi, rax	; pChild
	mov 	r8, r14		; j
	mov 	r9, r15		; i

LBL_COMPARE_WITH_valnow:
	cmp 	rdx, rbx
	jge		LBL_FOUND_LOCATION

; if here, current child is smaller, 
; so we push it up and go on
; vSieve[idxj_current] => rax
	mov     rax, qword ptr [rsi+r12*8]		
; vSieve[idxj_current][idxi_current] = pChild;
	mov     qword ptr [rax+r13*8], rdi	; 

; idx_current = idx_child
	mov		r12, r8
	mov		r13, r9

; new idx: 2*i + 1 
	shl		r8, 1
	shl		r9, 1
	inc		r9

; check for end of row
	cmp		r9, ROW_SIZE
	jge		LBL_ADVANCE_ON_NEXT_ROW

;  more probable branch here

LBL_TEST_FOR_END_LOOP:
; if idxj_child < last_idx_j continue
	cmp		r8, r10
	jl		LBL_START_LOOP
	jg		LBL_FOUND_LOCATION	; stop
; (here idxj are equal, we are in the last row)	
; if idxi_child < last_idx_i 
	cmp		r9, r11
	jg		LBL_FOUND_LOCATION
	jmp 	LBL_START_LOOP		; else continue

LBL_ADVANCE_ON_NEXT_ROW:
	inc		r8
	lea		r9, [r9-ROW_SIZE]	
	jmp		LBL_TEST_FOR_END_LOOP

LBL_FOUND_LOCATION:
; vSieve[idxj_current][idxi_current] = pRoot;
	;mov		rcx, pRoot
	mov		rcx, rbp
	mov     rax, qword ptr [rsi+r12*8]		; vSieve[idxj_current] 
	mov     qword ptr [rax+r13*8], rcx	; 
	
	pop rbp
	pop r15
	pop r14
	pop r13
	pop r12
	;pop r11
	;pop r10
	;pop r9
	;pop r8
	pop rsi
	pop rdi
	;pop rdx
	;pop rcx
	pop rbx
	;pop rax					
	ret

HeapifyASM endp

end