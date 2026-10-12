pro remchar,st,char	;Remove character
; ---------------------------------------------------------------
; Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

;+
; NAME:
;	REMCHAR
; PURPOSE:
;	Remove all appearances of character (char) from string (st)
;
; CALLING SEQUENCE:
;	REMCHAR, ST, CHAR
;
; INPUT-OUTPUT:
;	ST  - String from which character will be removed, scalar or vector
; INPUT:
;	CHAR- Single character to be removed from string or all elements of a
;		string array
;
; EXAMPLE:
;	If a = 'a,b,c,d,e,f,g' then
;
;	IDL> remchar,a, ','
;
;      will give a = 'abcdefg'
;
; REVISIONS HISTORY
;	Written D. Lindler October 1986
;	Test if empty string needs to be returned   W. Landsman  Feb 1991
;	Work on string arrays    W. Landsman   August 1997
;	Avoid 32 bit integer overflow K. Tolbert/W. Landsman Feb 2007
;-
 compile_opt idl2
 if N_params() LT 2 then begin
     print,'Syntax - REMCHAR, string, character'
     return
 endif

 bchar = byte(char) & bchar = bchar[0]          ;Convert character to byte

 for i = 0L,N_elements(st)-1 do  begin

 bst = byte(st[i])
 good = where( bst NE bchar, Ngood)
 if Ngood GT 0 then st[i] = string(bst[good]) else st[i] = ''

 endfor
 return
 end
