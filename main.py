# ============================================
# توابع کمکی
# ============================================

def print_matrix(mat, title="Matrix"):
    """چاپ ماتریس با فرمت مناسب"""
    print(f"\n{title}")
    for row in mat:
        print([round(x, 4) for x in row])

def copy_matrix(mat):
    """کپی عمیق از ماتریس"""
    return [row[:] for row in mat]

def swap_rows(mat, i, j):
    """جابجایی دو سطر"""
    if i != j:
        mat[i], mat[j] = mat[j], mat[i]

def scale_row(mat, i, scalar):
    """ضرب یک سطر در اسکالر"""
    mat[i] = [x * scalar for x in mat[i]]

def add_scaled_row(mat, dest, src, scalar):
    """اضافه کردن مضربی از سطر مبدأ به سطر مقصد"""
    for col in range(len(mat[0])):
        mat[dest][col] += scalar * mat[src][col]


# ============================================
# REF: تبدیل به فرم سطری پلکانی
# ============================================

def REF(matrix, verbose=True):
    """
    تبدیل ماتریس به فرم سطری پلکانی (Row Echelon Form)
    
    Parameters:
        matrix: ماتریس ورودی (لیست لیست‌ها)
        verbose: اگر True باشد، گام‌ها را چاپ می‌کند
    
    Returns:
        mat: ماتریس در فرم REF
        swaps: تعداد جابجایی‌های سطر
    """
    mat = copy_matrix(matrix)
    rows = len(mat)
    cols = len(mat[0]) if rows > 0 else 0
    
    row = 0
    swaps = 0
    
    if verbose:
        print_matrix(mat, "START: Initial matrix")
    
    for col in range(cols):
        if verbose:
            print(f"\n--- Checking column {col} ---")
        
        if row >= rows:
            if verbose:
                print("All rows processed, exiting loop")
            break
        
        # پیدا کردن سطر محور (pivot)
        pivot_row = None
        for r in range(row, rows):
            if abs(mat[r][col]) > 1e-10:
                pivot_row = r
                break
        
        if pivot_row is None:
            if verbose:
                print(f"Column {col}: No pivot found → moving to next column")
            continue
        
        if verbose:
            print(f"Pivot found at row {pivot_row}, column {col} (value = {mat[pivot_row][col]:.4f})")
        
        # جابجایی سطر محور به جایگاه درست
        if pivot_row != row:
            if verbose:
                print(f"Swapping row {row} with row {pivot_row}")
            swap_rows(mat, row, pivot_row)
            swaps += 1
            if verbose:
                print_matrix(mat, "After swap")
        
        # صفر کردن زیر محور
        for r in range(row + 1, rows):
            if abs(mat[r][col]) > 1e-10:
                factor = mat[r][col] / mat[row][col]
                if verbose:
                    print(f"Row {r} ← Row {r} - ({factor:.4f}) × Row {row}")
                add_scaled_row(mat, r, row, -factor)
                if verbose:
                    print_matrix(mat, f"After zeroing row {r}")
        
        row += 1
        if verbose:
            print(f"Moving to next row (row = {row})")
    
    if verbose:
        print("\n--- REF completed ---")
        print_matrix(mat, "REF matrix")
    
    return mat, swaps


# ============================================
# RREF: تبدیل به فرم سطری پلکانی تحویل‌یافته
# ============================================

def RREF(matrix, verbose=True):
    """
    تبدیل ماتریس به فرم سطری پلکانی تحویل‌یافته (Reduced Row Echelon Form)
    
    Parameters:
        matrix: ماتریس ورودی (لیست لیست‌ها)
        verbose: اگر True باشد، گام‌ها را چاپ می‌کند
    
    Returns:
        mat: ماتریس در فرم RREF
        swaps: تعداد جابجایی‌های سطر
    """
    # مرحله 1: تبدیل به REF
    mat, swaps = REF(matrix, verbose)
    
    rows = len(mat)
    cols = len(mat[0]) if rows > 0 else 0
    
    if verbose:
        print("\n" + "="*60)
        print("Starting RREF conversion (backward elimination)")
        print("="*60)
    
    # مرحله 2: از پایین به بالا پردازش کن
    for r in range(rows - 1, -1, -1):
        # پیدا کردن ستون پیشرو در سطر جاری
        pivot_col = None
        for c in range(cols):
            if abs(mat[r][c]) > 1e-10:
                pivot_col = c
                break
        
        if pivot_col is None:
            if verbose:
                print(f"\nRow {r} is all zeros → skipping")
            continue
        
        if verbose:
            print(f"\n--- Processing row {r} (pivot at column {pivot_col}) ---")
            print(f"Current pivot value = {mat[r][pivot_col]:.4f}")
        
        # نرمال‌سازی: پیشرو را به 1 تبدیل کن
        if abs(mat[r][pivot_col] - 1) > 1e-10:
            if verbose:
                print(f"Normalizing row {r}: divide by {mat[r][pivot_col]:.4f}")
            scale_row(mat, r, 1.0 / mat[r][pivot_col])
            if verbose:
                print_matrix(mat, f"After normalizing row {r}")
        
        # صفر کردن بالای پیشرو
        for r_up in range(r - 1, -1, -1):
            if abs(mat[r_up][pivot_col]) > 1e-10:
                factor = mat[r_up][pivot_col]
                if verbose:
                    print(f"Row {r_up} ← Row {r_up} - ({factor:.4f}) × Row {r}")
                add_scaled_row(mat, r_up, r, -factor)
                if verbose:
                    print_matrix(mat, f"After zeroing row {r_up}")
    
    if verbose:
        print("\n--- RREF completed ---")
        print_matrix(mat, "RREF matrix")
    
    return mat, swaps


# ============================================
# حل دستگاه معادلات خطی با روش گاوس-جردن
# ============================================

def solve_linear_system(A, b, verbose=True):
    """
    حل دستگاه معادلات خطی به فرم Ax = b با استفاده از روش گاوس-جردن
    
    Parameters:
        A: ماتریس ضرایب (لیست لیست‌ها)
        b: بردار سمت راست (لیست)
        verbose: اگر True باشد، گام‌ها را چاپ می‌کند
    
    Returns:
        x: جواب دستگاه (لیست)
    
    Raises:
        ValueError: اگر دستگاه جواب یکتا نداشته باشد
    """
    rows = len(A)
    cols = len(A[0])
    
    # بررسی سازگاری ابعاد
    if rows != len(b):
        raise ValueError("Dimensions of A and b do not match")
    
    # ساخت ماتریس افزوده
    augmented = [A[i][:] + [b[i]] for i in range(rows)]
    
    if verbose:
        print("\n" + "="*60)
        print("SOLVING LINEAR SYSTEM Ax = b USING GAUSS-JORDAN")
        print("="*60)
        print_matrix(augmented, "Augmented matrix [A|b]")
    
    # تبدیل به RREF
    rref_mat, _ = RREF(augmented, verbose)
    
    # استخراج جواب
    x = []
    
    # بررسی رتبه ماتریس و وجود جواب یکتا
    rank_A = 0
    rank_Aug = 0
    
    # محاسبه رتبه ماتریس ضرایب
    for i in range(rows):
        is_zero_row = True
        for j in range(cols):
            if abs(rref_mat[i][j]) > 1e-10:
                is_zero_row = False
                break
        if not is_zero_row:
            rank_A += 1
    
    # محاسبه رتبه ماتریس افزوده
    for i in range(rows):
        is_zero_row = True
        for j in range(cols + 1):
            if abs(rref_mat[i][j]) > 1e-10:
                is_zero_row = False
                break
        if not is_zero_row:
            rank_Aug += 1
    
    if verbose:
        print(f"\nRank(A) = {rank_A}")
        print(f"Rank([A|b]) = {rank_Aug}")
    
    if rank_A < rank_Aug:
        raise ValueError("System has no solution (inconsistent)")
    elif rank_A < cols:
        raise ValueError(f"System has infinite solutions (rank = {rank_A} < {cols})")
    
    # استخراج جواب (ماتریس باید به فرم [I|solution] باشد)
    for i in range(rows):
        # پیدا کردن ستون محور
        pivot_col = None
        for j in range(cols):
            if abs(rref_mat[i][j] - 1) < 1e-10:
                pivot_col = j
                break
        
        if pivot_col is not None:
            x.append(rref_mat[i][cols])
        else:
            x.append(0)  # برای سطرهای صفر
    
    if verbose:
        print("\n" + "="*60)
        print("SOLUTION:")
        for i, val in enumerate(x):
            print(f"x{i+1} = {val:.4f}")
    
    return x


# ============================================
# دترمینان با استفاده از REF
# ============================================

def determinant(matrix, verbose=True):
    """
    محاسبه دترمینان با استفاده از تبدیل به REF
    
    Parameters:
        matrix: ماتریس مربعی (لیست لیست‌ها)
        verbose: اگر True باشد، گام‌ها را چاپ می‌کند
    
    Returns:
        det: مقدار دترمینان
    """
    mat = copy_matrix(matrix)
    rows = len(mat)
    cols = len(mat[0]) if rows > 0 else 0
    
    # بررسی مربعی بودن ماتریس
    if rows != cols:
        raise ValueError("Determinant only defined for square matrices")
    
    if verbose:
        print_matrix(mat, "START: Original matrix for determinant")
    
    det = 1.0
    swaps = 0
    row = 0
    
    for col in range(cols):
        if row >= rows:
            break
        
        # پیدا کردن سطر محور
        pivot_row = None
        for r in range(row, rows):
            if abs(mat[r][col]) > 1e-10:
                pivot_row = r
                break
        
        # اگر هیچ محوری پیدا نشد، دترمینان صفر است
        if pivot_row is None:
            if verbose:
                print(f"\nColumn {col}: No pivot found → determinant = 0")
            return 0.0
        
        # جابجایی در صورت نیاز
        if pivot_row != row:
            if verbose:
                print(f"\nSwapping row {row} with row {pivot_row}")
            swap_rows(mat, row, pivot_row)
            swaps += 1
            if verbose:
                print_matrix(mat, "After swap")
        
        # صفر کردن زیر محور (بدون تغییر در دترمینان)
        for r in range(row + 1, rows):
            if abs(mat[r][col]) > 1e-10:
                factor = mat[r][col] / mat[row][col]
                if verbose:
                    print(f"Row {r} ← Row {r} - ({factor:.4f}) × Row {row}")
                add_scaled_row(mat, r, row, -factor)
                if verbose:
                    print_matrix(mat, f"After elimination")
        
        row += 1
    
    # محاسبه دترمینان از قطر اصلی
    if verbose:
        print("\n--- REF completed for determinant ---")
        print_matrix(mat, "REF matrix")
    
    for i in range(rows):
        det *= mat[i][i]
        if verbose:
            print(f"Diagonal[{i}] = {mat[i][i]:.4f}, product so far = {det:.4f}")
    
    if swaps % 2 == 1:
        det = -det
        if verbose:
            print(f"\nSwaps count = {swaps} (odd) → multiplying by -1")
    
    if verbose:
        print(f"\nFinal determinant = {det:.4f}")
    
    return det


# ============================================
# محاسبه ماتریس وارون
# ============================================

def inverse_matrix(matrix, verbose=True):
    """
    محاسبه ماتریس وارون با استفاده از روش گاوس-جردن
    
    Parameters:
        matrix: ماتریس مربعی (لیست لیست‌ها)
        verbose: اگر True باشد، گام‌ها را چاپ می‌کند
    
    Returns:
        inv: ماتریس وارون
    
    Raises:
        ValueError: اگر ماتریس وارون‌پذیر نباشد
    """
    rows = len(matrix)
    cols = len(matrix[0])
    
    # بررسی مربعی بودن ماتریس
    if rows != cols:
        raise ValueError("Inverse only defined for square matrices")
    
    # محاسبه دترمینان
    det = determinant(matrix, verbose=False)
    
    if abs(det) < 1e-10:
        raise ValueError(f"Matrix is singular (determinant = {det:.6f}), cannot compute inverse")
    
    if verbose:
        print("\n" + "="*60)
        print(f"COMPUTING INVERSE MATRIX (det = {det:.4f})")
        print("="*60)
        print_matrix(matrix, "Original matrix")
    
    # ساخت ماتریس افزوده [A|I]
    augmented = []
    for i in range(rows):
        row = matrix[i][:] + [1.0 if j == i else 0.0 for j in range(cols)]
        augmented.append(row)
    
    if verbose:
        print_matrix(augmented, "Augmented matrix [A|I]")
    
    # اعمال روش گاوس-جردن
    rref_mat, _ = RREF(augmented, verbose)
    
    # استخراج ماتریس وارون (نیمه راست ماتریس افزوده)
    inv = []
    for i in range(rows):
        inv_row = []
        for j in range(cols):
            inv_row.append(rref_mat[i][cols + j])
        inv.append(inv_row)
    
    if verbose:
        print("\n" + "="*60)
        print("INVERSE MATRIX:")
        print_matrix(inv, "Inverse")
        
        # تست صحت: A × A^(-1) باید برابر I باشد
        print("\n" + "="*60)
        print("VERIFICATION: A × A^(-1)")
        result = matrix_multiply(matrix, inv)
        print_matrix(result, "A × A^(-1)")
    
    return inv


# ============================================
# توابع کمکی اضافی
# ============================================

def matrix_multiply(A, B):
    """ضرب دو ماتریس"""
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])
    
    if cols_A != rows_B:
        raise ValueError("Cannot multiply matrices: incompatible dimensions")
    
    result = [[0.0 for _ in range(cols_B)] for _ in range(rows_A)]
    
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                result[i][j] += A[i][k] * B[k][j]
    
    return result


def get_matrix_from_user():
    """دریافت ماتریس از کاربر"""
    print("\n" + "="*60)
    print("ENTER MATRIX")
    print("="*60)
    
    while True:
        try:
            rows = int(input("Number of rows: "))
            cols = int(input("Number of columns: "))
            
            if rows <= 0 or cols <= 0:
                print("Rows and columns must be positive!")
                continue
            
            matrix = []
            print(f"\nEnter matrix elements row by row:")
            
            for i in range(rows):
                while True:
                    try:
                        row_str = input(f"Row {i+1} (enter {cols} numbers separated by space): ")
                        row = [float(x) for x in row_str.split()]
                        
                        if len(row) != cols:
                            print(f"Please enter exactly {cols} numbers!")
                            continue
                        
                        matrix.append(row)
                        break
                    except ValueError:
                        print("Invalid input! Please enter numbers only.")
            
            return matrix
        
        except ValueError:
            print("Invalid number format!")


def get_vector_from_user(n):
    """دریافت بردار از کاربر"""
    print(f"\nEnter vector b (size {n}):")
    
    while True:
        try:
            vec_str = input(f"Enter {n} numbers separated by space: ")
            vec = [float(x) for x in vec_str.split()]
            
            if len(vec) != n:
                print(f"Please enter exactly {n} numbers!")
                continue
            
            return vec
        except ValueError:
            print("Invalid input! Please enter numbers only.")


# ============================================
# منوی اصلی برنامه
# ============================================

def main():
    """منوی اصلی برنامه"""
    while True:
        print("\n" + "="*60)
        print("LINEAR ALGEBRA PROJECT - GAUSS-JORDAN ELIMINATION")
        print("="*60)
        print("1. Calculate Determinant")
        print("2. Solve Linear System (Ax = b)")
        print("3. Compute Inverse Matrix")
        print("4. Exit")
        print("="*60)
        
        choice = input("Choose an option (1-4): ")
        
        if choice == '1':
            print("\n" + "-"*60)
            print("DETERMINANT CALCULATION")
            print("-"*60)
            
            matrix = get_matrix_from_user()
            
            if len(matrix) != len(matrix[0]):
                print(f"\nError: Matrix is not square ({len(matrix)}×{len(matrix[0])})")
                print("Determinant only defined for square matrices!")
            else:
                verbose = input("\nShow steps? (y/n): ").lower() == 'y'
                
                try:
                    det = determinant(matrix, verbose)
                    print(f"\n{'='*60}")
                    print(f"DETERMINANT = {det:.6f}")
                    print(f"{'='*60}")
                except Exception as e:
                    print(f"\nError: {e}")
        
        elif choice == '2':
            print("\n" + "-"*60)
            print("SOLVE LINEAR SYSTEM Ax = b")
            print("-"*60)
            
            A = get_matrix_from_user()
            n = len(A)
            
            # بررسی مربعی بودن ماتریس
            if len(A[0]) != n:
                print(f"\nError: Matrix must be square for solving system!")
                continue
            
            b = get_vector_from_user(n)
            
            verbose = input("\nShow steps? (y/n): ").lower() == 'y'
            
            try:
                x = solve_linear_system(A, b, verbose)
                print(f"\n{'='*60}")
                print("SOLUTION:")
                for i, val in enumerate(x):
                    print(f"x{i+1} = {val:.6f}")
                print(f"{'='*60}")
            except Exception as e:
                print(f"\nError: {e}")
        
        elif choice == '3':
            print("\n" + "-"*60)
            print("INVERSE MATRIX CALCULATION")
            print("-"*60)
            
            matrix = get_matrix_from_user()
            
            if len(matrix) != len(matrix[0]):
                print(f"\nError: Matrix is not square ({len(matrix)}×{len(matrix[0])})")
                print("Inverse only defined for square matrices!")
            else:
                verbose = input("\nShow steps? (y/n): ").lower() == 'y'
                
                try:
                    inv = inverse_matrix(matrix, verbose)
                    print(f"\n{'='*60}")
                    print("INVERSE MATRIX:")
                    for row in inv:
                        print([round(x, 6) for x in row])
                    print(f"{'='*60}")
                except Exception as e:
                    print(f"\nError: {e}")
        
        elif choice == '4':
            print("\nExiting program. Goodbye!")
            break
        
        else:
            print("\nInvalid choice! Please select 1-4.")


# اجرای برنامه
if __name__ == "__main__":
    main()
