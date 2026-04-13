from itertools import permutations


def prepare_matrix(matrixA, vectorB):
    n = len(matrixA)

    if any(len(row) != n for row in matrixA):
        raise ValueError("המטריצה חייבת להיות ריבועית")

    if len(vectorB) != n:
        raise ValueError("גודל וקטור B חייב להתאים לגודל המטריצה")

    A = [[float(x) for x in row] for row in matrixA]

    b = []
    for item in vectorB:
        if isinstance(item, list):
            b.append(float(item[0]))
        else:
            b.append(float(item))

    def is_diagonally_dominant(mat):
        for i in range(len(mat)):
            diagonal = abs(mat[i][i])
            others_sum = sum(abs(mat[i][j]) for j in range(len(mat)) if j != i)
            if diagonal < others_sum:
                return False
        return True

    if is_diagonally_dominant(A):
        return A, b, True, True

    for perm in permutations(range(n)):
        new_A = [A[i][:] for i in perm]
        new_b = [b[i] for i in perm]

        if is_diagonally_dominant(new_A):
            return new_A, new_b, True, False

    return A, b, False, False


def jacobi(matrixA, vectorB, tolerance=0.00001, max_iterations=100):
    n = len(matrixA)
    x_old = [0.0] * n
    x_new = [0.0] * n

    for iteration in range(1, max_iterations + 1):
        for i in range(n):
            if matrixA[i][i] == 0:
                raise ZeroDivisionError("אי אפשר לחלק ב-0 על האלכסון הראשי")

            sigma = 0.0
            for j in range(n):
                if j != i:
                    sigma += matrixA[i][j] * x_old[j]

            x_new[i] = (vectorB[i] - sigma) / matrixA[i][i]

        print(f"איטרציה {iteration}: {x_new}")

        max_error = 0.0
        for i in range(n):
            current_error = abs(x_new[i] - x_old[i])
            if current_error > max_error:
                max_error = current_error

        if max_error < tolerance:
            return x_new[:], iteration, True

        x_old = x_new[:]

    return x_new[:], max_iterations, False


def gauss_seidel(matrixA, vectorB, tolerance=0.00001, max_iterations=100):
    n = len(matrixA)
    x = [0.0] * n

    for iteration in range(1, max_iterations + 1):
        old_x = x[:]

        for i in range(n):
            if matrixA[i][i] == 0:
                raise ZeroDivisionError("אי אפשר לחלק ב-0 על האלכסון הראשי")

            sum_before = 0.0
            sum_after = 0.0

            for j in range(i):
                sum_before += matrixA[i][j] * x[j]

            for j in range(i + 1, n):
                sum_after += matrixA[i][j] * old_x[j]

            x[i] = (vectorB[i] - sum_before - sum_after) / matrixA[i][i]

        print(f"איטרציה {iteration}: {x}")

        max_error = 0.0
        for i in range(n):
            current_error = abs(x[i] - old_x[i])
            if current_error > max_error:
                max_error = current_error

        if max_error < tolerance:
            return x[:], iteration, True

    return x[:], max_iterations, False


def run_method(method_name, method_function, A, B, has_dominant_diagonal):
    print(f"\n--- {method_name} ---")

    solution, iterations, converged = method_function(A, B)

    print("פתרון סופי:", solution)
    print("מספר איטרציות:", iterations)

    if has_dominant_diagonal:
        if converged:
            print("השיטה התכנסה.")
        else:
            print("השיטה לא התכנסה במסגרת מספר האיטרציות.")
    else:
        if converged:
            print("למרות שאין אלכסון דומיננטי, התוצאות הן:", solution)
        else:
            print("המערכת אינה מתכנסת")


# ======================
# תוכנית ראשית
# ======================

matrixA = [[4, 2, 0],
           [2, 10, 4],
           [0, 4, 5]]

vectorB = [[2],
           [6],
           [5]]

A, B, has_dominant_diagonal, was_originally_dominant = prepare_matrix(matrixA, vectorB)

print("המטריצה שבה נשתמש:")
for row in A:
    print(row)

print("וקטור B:")
print(B)

if has_dominant_diagonal:
    if was_originally_dominant:
        print("המטריצה הייתה בעלת אלכסון דומיננטי מההתחלה.")
    else:
        print("בוצע סידור שורות כדי לקבל אלכסון דומיננטי.")
else:
    print("לא ניתן להביא את המטריצה לצורה עם אלכסון דומיננטי.")
    print("נריץ בכל זאת מספר איטרציות ונבדוק האם יש התכנסות.")

run_method("שיטת יעקובי", jacobi, A, B, has_dominant_diagonal)
run_method("שיטת גאוס-זיידל", gauss_seidel, A, B, has_dominant_diagonal)
