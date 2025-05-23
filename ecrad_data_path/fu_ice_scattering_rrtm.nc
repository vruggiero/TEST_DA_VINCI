CDF       
      band_lw       band_sw       coeff_sw   
   coeff_lw            comment      This file provides a parameterization of ice particle scattering in the longwave and shortwave RRTM bands,
using the formulation of Fu (J. Climate, 1996, 9, 2058-2082) in the shortwave and
Fu et al. (J. Climate, 1998, 11, 2223-2237) in the longwave. If De is Fu's ice effective
size (microns), and the 1-based arrays of coefficients are denoted p, then in the shortwave:
  mass extinction coefficient (m2/g) = p[1] + p[2]/De,
  single scattering albedo = 1 - (p[3] + p[4]*De + p[5]*De^2 + p[6]*De^3), and
  asymmetry factor = p[7] + p[8]*De + p[9]*De^2 + p[10]*De^3,
and in the longwave:
  mass extinction coefficient (m2/g) = p[1] + p[2]/De + p[3]/De^2,
  single scattering albedo = 1 - (p[4] + p[5]*De + p[6]*De^2 + p[7]*De^3), and
  asymmetry factor = p[8] + p[9]*De + p[10]*De^2 + p[11]*De^3,
where here mass extinction coefficient is the total extinction cross section per unit mass of ice.          wavenumber1_lw                  	long_name         (Lower bound wavenumber for longwave band   units         cm-1      @  �   wavenumber2_lw                  	long_name         (Upper bound wavenumber for longwave band   units         cm-1      @  0   wavenumber1_sw                 	long_name         )Lower bound wavenumber for shortwave band      units         cm-1      8  p   wavenumber2_sw                 	long_name         )Upper bound wavenumber for shortwave band      units         cm-1      8  �   coeff_lw                   	long_name         $Longwave ice scattering coefficients     �  �   coeff_sw                  	long_name         %Shortwave ice scattering coefficients        0  
�A   C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� EK  E"� EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp DM  EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp GCP E"� ;�5O@���^��?c	<��B��ӵ5���>�g&<BW�9�4�'
;Y��@&|���Ħ>M_<��(�7�r5�]?0k�;�螸�w@4��]�Z�@�f�����>��<���J�?50�??9�O;.��Q��2��b�H@e���+�?t�Q<V���<5DQ??<�6;�:Ƹ9_4#�(��K�@C_���]?�����|�dq�4 *�?Ld;�Rz�;��45���ݹ�@1%�K�?��;[��\��4oΙ?]�\;7]���E4�c�	-�@]�r��?
{�<���W�5E��?d 2:�w��k�3{y0���@Y� �P?a�<����V�p5I��?\g�;5W���E3^U���I@E���x��?7%<��:�C�Y5>�?Z/ ;%co��I�3����F߰@9�N����?b��<H�����5!�?_�`; ����G3����F߰@9�N����?b��<H�����5!�?_�`; ����G3����!��@5R��s�?>.j<Q�[���5g?X�;'�ⷷD3���P�@4�?�>6>ފ5<���*5+5J?K�u;D�+����3f}ɻP�@4�?�>6>ފ5<���*5+5J?K�u;D�+����3f}ɻP�@4�?�>6>ފ5<���*5+5J?K�u;D�+����3f}ɻP�@4�?�>6>ފ5<���*5+5J?K�u;D�+����3f}�9D��@ �>I�;�~6�]��4=�c?BY�;���Kp429��B@J�>���7��m��)u1��?kj�:�p���t2�I: P?@7v:��:��ضd-r2��?F�n:�0��#2����z�@!�A;%��:�zn�bn2Ʊc?B=O;���DW�2���z�@!�A;%��:�zn�bn2Ʊc?B=O;���DW�2���S
@$�/�:�oF���2��?@��; I�+�D2�Fa8���@!���;8�<Ѳ�>q.0�~?@S�:�)���v71���8��T@!	B�)��5�}�.�bު՘�?@��:��Ŷ��}9)��@ z:��4������+��??��:�5ζA�v:	��F�@!Q�5|62ꮢ.Zb�Qk#??��:i�����2<�@���U@"q�U+<3��v���R+�J?>Y:n4����5�-���@"������4,E��b �+z?<�:kG6�i��Gp�w��@"qa����4d�����O+�j?;��:p�6���9�Y����@!ܜ>�3	;��͸j\f4S��?L��;��n�WB�4F�