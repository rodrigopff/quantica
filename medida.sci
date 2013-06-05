//***************************************************************
//Calcula a formula de medida de um q-bit num estado emaranhado
//***************************************************************

//######################
//SEÇÃO DE CONSTANTES  #
//######################

//Matrizes de Pauli
omegaX= [0 1;1 0]
omegaY= [0 -%i;%i 0]
omegaZ= [1 0;0 -1]

//Mapa de Matrizes Pauli indexada pela posição de 1 a 3 
listaMatrizesPauli=list(omegaX,omegaY,omegaZ)

//###############################
// SEÇÃO DE FUNÇÕES DE CÁLCULO  #
//###############################

//Função que extrai a raiz quadrada de uma matriz de um estado quântico emaranhado 
//Criar um tratamento para tratar a diagonalização de matrizes que é pré-requisito para que uma Matriz possua raiz quandrada.
function [raizMatriz] = raizRO (roab)
    raizMatriz=sqrtm(roab)
endfunction

//Função que calcula o produto tensorial entre os estados sigma_xyz_A e a matriz identidade Ib 

//function [produto_Tensorial] = calculaTensorial(sigmaXYZ_A, Ib)
function [produtoTensorial] = calculaTensorial(listaMatrizesTensorial)
    //produto_Tensorial
    produtoTensorial=1 //elemento neutro da multiplicação como escalar
    for i= 1:length(listaMatrizesTensorial) 
          produtoTensorial=produtoTensorial.*.listaMatrizesTensorial(i)
    end
endfunction

//Função que aplica sigma-Pauli-XYZ ao primeiro par do estado quântico 

function [sigmaXYZ_A] = aplicaPauliA (pauliMatriz,estadoQuantico)
    sigmaXYZ_A = pauliMatriz*estadoQuantico    
endfunction
  
//Função que multiplica N matrizes 
function [produtoMatriz]= produtoMatrizes (listaMatrizesProduto)
    produtoMatriz=1 //elemento neutro da multiplicação como escalar
    for i= 1:length(listaMatrizesProduto) 
          produtoMatriz=produtoMatriz*listaMatrizesProduto(i)
    end
endfunction

//Função que calcula a matriz Identidade de um vetor de estado quantico
function [Id] = Ident(vetorEstado)
    Id=eye(length(vetorEstado), length(vetorEstado))
endfunction

//Função que itera sobre as aplicações de Pauli com x=1,y=2 e z=3
//function () = 
//endfunction

//Função que calcula a matriz Wab

// parametro de entrada:
// matriz |phi><phi| , estado quantico do primeiro q-bit, estado quantico do segundo q-bit
// parametro

//Caucula a matrizW

function [listaResultado] = calculaMatrizW (roab,estadoQuanticoA,estadoQuanticoB)
    //Indexar os indices i,j correspondentes às operações x,y,z
    //Serão utilizado dois dois for's aninhados. Talvez isto possa ser melhorado realizando uma refatoração de código
    listaResultados=[]
    //Itera sobre i variando de x a z (de 1 a 3). Verifica se esta variação pode ser maior que 3 para outras medidas para parametrizar os laços abaixo    
    for i=1:3 
        for j=1:3                                 
   matrizW=trace(produtoMatrizes(list(raizRO(roab),calculaTensorial(list(listaMatrizesPauli(i),Ident(estadoQuanticoB))),raizRO(roab),calculaTensorial(list(listaMatrizesPauli(j),Ident(estadoQuanticoB))))))
   listaResultado(i,j)=matrizW    

    //matrizW=trace(produtoMatrizes(list(raizRO(roab),calculaTensorial(list(aplicaPauliA(listaMatrizesPauli(i),estadoQuanticoA),Ident(estadoQuanticoB))),raizRO(roab),calculaTensorial(list(aplicaPauliA(listaMatrizesPauli(j),estadoQuanticoA),Ident(estadoQuanticoB))))))
        end        
    end        
    
    //trace (ProdutoMatrizes(raizRO(roab),calculaTensorial(,Ident(estadoQuanticoB)),raizRO(roab),calculaTensorial(,Ident(estadoQuanticoB))))
    
    //trace(raizRO(roab)*calculaTensorial(,Ident(estadoQuanticoB))*raizRO(roab)*calculaTensorial(,Ident(estadoQuanticoB)))    
endfunction



//function[tes1] = teste1 ()
//    tes1=teste2('Chamando de teste1 ......')
//endfunction

//function[tes2]  = teste2 (b)
//    tes2=b
//endfunction    



//#######################################
// SEÇÃO DE FUNÇÕES DE ENTRADA E SAÍDA  #
//#######################################

