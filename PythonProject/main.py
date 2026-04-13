from itertools import permutations


def prepare_matrix(matrixA, vectorB):
    n = len(matrixA)

    if any(len(row) != n for row in matrixA):
        raise ValueError("המטריצה חייבת להיות ריבועית")

    if len(vectorB) != n:
        raise ValueError("גודל וקטור B חייב להתאים לגודל המטריצה")

    b = []
    for item in vectorB:
        if isinstance(item, list):
            b.append(float(item[0]))
        else:
            b.append(float(item))

    A = []
    for row in matrixA:
        A.append([float(x) for x in row])

    def is_diagonally_dominant(mat):
        for i in range(len(mat)):
            diagonal = abs(mat[i][i])
            others_sum = 0
            for j in range(len(mat)):
                if i != j:
                    others_sum += abs(mat[i][j])
            if diagonal < others_sum:
                return False
        return True


    if is_diagonally_dominant(A):
        return A, b, True, True


    for perm in permutations(range(n)):
        new_A = []
        new_b = []
        for i in perm:
            new_A.append(A[i][:])
            new_b.append(b[i])

        if is_diagonally_dominant(new_A):
            return new_A, new_b, True, False

    return A, b, False, False


def gauss_seidel(matrixA, vectorB, tolerance=0.00001, max_iterations=100):
    n = len(matrixA)
    x = [0.0] * n

    for iteration in range(1, max_iterations + 1):
        old_x = x[:]

        for i in range(n):
            sum_before = 0.0
            sum_after = 0.0

            for j in range(i):
                sum_before += matrixA[i][j] * x[j]

            for j in range(i + 1, n):
                sum_after += matrixA[i][j] * old_x[j]

            if matrixA[i][i] == 0:
                raise ZeroDivisionError("אי אפשר לחלק ב-0 על האלכסון הראשי")

            x[i] = (vectorB[i] - sum_before - sum_after) / matrixA[i][i]

        print(f"איטרציה {iteration}: {x}")

        max_error = 0.0
        for i in range(n):
            current_error = abs(x[i] - old_x[i])
            if current_error > max_error:
                max_error = current_error

        if max_error < tolerance:
            return x, iteration, True

    return x, max_iterations, False


# תוכנית ראשית
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
print()

if has_dominant_diagonal:
    if was_originally_dominant:
        print("המטריצה היא עם אלכסון דומיננטי.")
    else:
        print("המטריצה לא הייתה עם אלכסון דומיננטי, בוצעה החלפת שורות למטריצה עם אלכסון דומיננטי.")

    result, iterations, converged = gauss_seidel(A, B)

    if converged:
        print()
        print("התוצאות הן:")
        for i in range(len(result)):
            print(f"x{i + 1} = {result[i]}")
        print(f"מספר האיטרציות שהתבצעו: {iterations}")
    else:
        print()
        print("המערכת אינה מתכנסת.")
        print(f"מספר האיטרציות שהתבצעו: {iterations}")

else:
    print("לא התקבלה מטריצה עם אלכסון דומיננטי גם לאחר החלפת שורות.")
    print("ננסה בכל זאת לפתור בעזרת גאוס זיידל עד למספר איטרציות מקסימלי.")
    print()

    result, iterations, converged = gauss_seidel(A, B)

    print()
    if converged:
        print("למטריצה שאין אלכסון דומיננטי התוצאות הן:")
        for i in range(len(result)):
            print(f"x{i + 1} = {result[i]}")
        print(f"מספר האיטרציות שהתבצעו: {iterations}")
    else:
        print("המערכת אינה מתכנסת.")
        print(f"מספר האיטרציות שהתבצעו: {iterations}")
