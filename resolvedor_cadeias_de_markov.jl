using LinearAlgebra

# Função para verificar se a matriz é estocástica por colunas
function verificar_estocastica_por_colunas(T::Matrix{Float64})
    n_rows, n_cols = size(T)

    if n_rows != n_cols
        println("Erro: A matriz de transição deve ser quadrada.")
        return false
    end

    for j in 1:n_cols

        soma_coluna = sum(T[:, j])
        if abs(soma_coluna - 1.0) > 1e-5 

            println("Erro: A soma da coluna $j é $soma_coluna, mas deveria ser 1.")
            return false

        end

    end

    return true

end

# Função para obter a matriz de transição inserida manualmente
function obter_matriz_manual()::Union{Matrix{Float64}, Nothing}

    n = 0

    while n <= 0

        try

            println("Digite o número de estados (tamanho da matriz, ex: 2, 3, 4...):")
            n = parse(Int, readline())

            if n <= 0

                println("O número de estados deve ser positivo.")

            end

        catch

            println("Entrada inválida. Digite um número inteiro positivo.")
            n = 0 # Reseta para continuar o loop

        end

    end
    
    T_temp = Matrix{Float64}(undef, n, n)
    println("\nInsira os valores da matriz de transição linha por linha.")
    println("\nOs valores em cada linha devem ser separados por espaço.")
    println("\nLembre-se: a SOMA DE CADA COLUNA deve ser 1.")

    while true # Loop para tentar inserir a matriz até ser válida

        valida = true

        for i in 1:n

            while true # Loop para tentar ler uma linha válida

                try

                    println("Linha $i (insira $n valores separados por espaço):")
                    entrada = split(readline())

                    if length(entrada) != n

                        println("Erro: Digite exatamente $n valores.")
                        continue # Tenta ler a linha novamente

                    end

                    valores_linha = parse.(Float64, entrada) # Verifica se há valores negativos

                    if any(x -> x < 0, valores_linha)
                        println("Erro: As probabilidades não podem ser negativas.")
                        continue # Tenta ler a linha novamente
                    end

                    T_temp[i, :] = valores_linha
                    break # Sai do loop de leitura da linha

                catch e 
                    println("Erro ao ler os valores da linha $i. Verifique se são números válidos separados por espaço. Detalhe: $e")
                    # Continua no loop para tentar ler a linha novamente

                end

            end # Fim do loop de leitura da linha
        
        end # Fim do loop pelas linhas

        # Verifica se a matriz é estocástica por colunas APÓS ler todas as linhas
        if verificar_estocastica_por_colunas(T_temp)
            return T_temp # Matriz válida
        else
            
            println("Erro: A matriz inserida não é estocástica por colunas (a soma de cada coluna deve ser 1). Por favor, insira a matriz novamente.")
            # Continua no loop para tentar inserir a matriz novamente
        end
    end # Fim do loop de tentativa de inserção
end

# Função para calcular a distribuição estacionária (CORRIGIDA)
function calcular_distribuicao_estacionaria(T::Matrix{Float64})
    n = size(T, 1)
    try
        # Montar o sistema A * π = b
        # A partir de (T - I)π = 0 e sum(π) = 1
        A_temp = T - I 
        A = Matrix{Float64}(undef, n, n)
        
        # Copiar as primeiras n-1 linhas de (T - I) para A
        # Estas representam as equações (T-I)π = 0 (que são linearmente dependentes)
        if n > 1
            A[1:n-1, :] = A_temp[1:n-1, :]
        end
        
        # Substituir a última linha de A pela equação de restrição sum(π) = 1
        A[n, :] = ones(1, n)
        
        # Montar o vetor b
        # As primeiras n-1 equações são = 0
        # A última equação (restrição) é = 1
        b = zeros(n)
        b[n] = 1.0

        # Resolver sistema linear A * π = b
        π = A \ b

        # Verifica se a solução é válida (soma 1 e não negativa)
        if abs(sum(π) - 1.0) > 1e-5 || any(x -> x < -1e-6, π)
             println("\n  Atenção: A solução numérica encontrada pode não ser uma distribuição de probabilidade válida.")
             println("  π calculado: $(round.(π, digits=6))")
             println("  Tentando método alternativo (Autovetor)...")
             
             # Método alternativo: Autovetor associado ao autovalor 1
             try 
                vals, vecs = eigen(T)
                idx = findfirst(v -> isapprox(v, 1.0), vals)
                if idx !== nothing
                    π_eigen = real.(vecs[:, idx])
                    π_eigen ./= sum(π_eigen) # Normaliza para somar 1
                    if abs(sum(π_eigen) - 1.0) < 1e-5 && !any(x -> x < -1e-6, π_eigen)
                        π = π_eigen # Usa a solução do autovetor se for válida
                        println("  Solução encontrada via autovetor.")
                    else
                        println("  Solução via autovetor também inválida.")
                    end
                else
                    println("  Não foi encontrado autovalor próximo de 1.")
                end
             catch e_eigen
                println("  Erro ao calcular autovetores: $e_eigen")
             end
        end

        println("\n  Distribuição estacionária (π) encontrada:")
        for i in 1:n
            println("    π[$i] = $(round(π[i], digits=6))")
        end

        println("\n  Verificação (T * π):")
        verificacao = T * π
        println("    $(round.(verificacao, digits=6))")
        println("    (Deveria ser aproximadamente igual a π)")

    catch e
        println("  Ocorreu um erro ao calcular a distribuição estacionária: ", e)
        if isa(e, LinearAlgebra.SingularException)
            println("  A matriz do sistema linear pode ser singular. Verifique se a cadeia de Markov é regular/ergódica ou se há estados absorventes.")
        end
    end
end

# Função para calcular o estado após n passos
function calcular_estado_n_passos(T::Matrix{Float64})
    n_estados = size(T, 1)
    p1 = nothing
    
    # Obter vetor de estado inicial p^(1)
    while p1 === nothing
        println("\n  Digite o vetor de estado inicial p^(1) (com $n_estados valores separados por espaço):")
        println("  A soma dos valores deve ser 1.")
        try
            entrada = split(readline())
            if length(entrada) != n_estados
                println("  Erro: Digite exatamente $n_estados valores.")
                continue
            end
            p1_temp = parse.(Float64, entrada)
            if abs(sum(p1_temp) - 1.0) > 1e-5
                 println("  Erro: A soma das probabilidades no vetor inicial deve ser 1 (Soma atual: $(sum(p1_temp))).")
                 continue
            end
             if any(x -> x < 0, p1_temp)
                println("  Erro: As probabilidades não podem ser negativas.")
                continue
            end
            p1 = p1_temp # Vetor válido
        catch e
            println("  Erro ao ler o vetor inicial. Verifique se são números válidos separados por espaço. Detalhe: $e")
        end
    end

    # Obter número de passos n
    num_passos = 0
    while num_passos <= 0
        println("\n  Digite o número de passos (n) para a previsão (deve ser um inteiro positivo):")
        try
            num_passos = parse(Int, readline())
            if num_passos <= 0
                println("  O número de passos deve ser positivo.")
            end
        catch e
            println("  Entrada inválida. Digite um número inteiro positivo. Detalhe: $e")
            num_passos = 0 # Reseta para continuar o loop
        end
    end

    # Calcular T^n
    try
        Tn = T^num_passos
        println("\n  Matriz T^$num_passos:")
        # Imprimir a matriz formatada
        for i = 1:n_estados
            print("    ")
            for j = 1:n_estados
                print("$(round(Tn[i,j], digits=6)) ")
            end
            println()
        end

        # Calcular p^(n) = T^n * p^(1)
        pn = Tn * p1

        println("\n  Vetor de estado após $num_passos passos (p^($num_passos) = T^$num_passos * p^(1)):")
        for i in 1:n_estados
            println("    p^($num_passos)[$i] = $(round(pn[i], digits=6))")
        end
         println("  Soma: $(round(sum(pn), digits=6))") # Deve ser 1

    catch e
        println("  Ocorreu um erro ao calcular T^n ou p^(n): ", e)
    end
end

# Função para executar os cálculos para uma dada matriz
function run_calculations_for_matrix(T::Matrix{Float64}, titulo_cenario::String, show_continue_option::Bool)
    println("\n=================================================")
    println(titulo_cenario)
    println("-------------------------------------------------")
    println("Matriz de Transição (T):")
    n_estados = size(T,1)
    for i = 1:n_estados
        print("  ")
        for j = 1:n_estados
            print("$(T[i,j]) ")
        end
        println()
    end
    println("-------------------------------------------------")

    while true
        println("\nO que você deseja calcular para este cenário?")
        println("  1: Distribuição Estacionária (longo prazo)")
        println("  2: Estado após n passos")
        if show_continue_option
            println("  3: Continuar para o próximo cenário")
        else
            println("  3: Sair")
        end
        
        try
            escolha_calculo = parse(Int, readline())

            if escolha_calculo == 1
                calcular_distribuicao_estacionaria(T)
            elseif escolha_calculo == 2
                calcular_estado_n_passos(T)
            elseif escolha_calculo == 3
                return true # Indica para continuar (ou sair se for o último)
            else
                println("Opção inválida. Escolha 1, 2 ou 3.")
            end
        catch e
             println("Entrada inválida. Digite 1, 2 ou 3. Detalhe: $e")
        end
        println("-------------------------------------------------") # Separador para clareza
    end
end

# Função principal
function main()
    println("=================================================")
    println("       Resolvedor de Cadeias de Markov")
    println("=================================================")

    modo_manual = nothing
    while modo_manual === nothing
        println("Deseja inserir manualmente a matriz de transição? (s/n): ")
        resp = lowercase(strip(readline()))
        if resp == "s"
            modo_manual = true
        elseif resp == "n"
            modo_manual = false
        else
            println("Resposta inválida. Digite 's' ou 'n'.")
        end
    end

    if modo_manual
        println("\n--- Modo de Entrada Manual ---")
        T_manual = obter_matriz_manual()
        if T_manual !== nothing
            run_calculations_for_matrix(T_manual, "Cenário: Matriz Inserida Manualmente", false) # false = Sair no final
        else
            println("Erro ao obter matriz manual. Encerrando.")
        end
    else
        println("\n--- Modo de Comparação Automática ---")
        # Define as matrizes para os cenários
        T_naturais = [0.6 0.1;
                      0.4 0.9]
        T_remedios = [0.3 0.05;
                      0.7 0.95]
        
        println("Serão analisados dois cenários pré-definidos:")
        println("1. Cuidados Naturais")
        println("2. Uso de Remédios")

        # Executa para o primeiro cenário
        if run_calculations_for_matrix(T_naturais, "Cenário 1: Cuidados Naturais", true) # true = Continuar
            # Executa para o segundo cenário
            run_calculations_for_matrix(T_remedios, "Cenário 2: Uso de Remédios", false) # false = Sair
        end
    end

    println("\n=================================================")
    println("Programa encerrado.")
    println("=================================================")
end

# Executa o programa
main()