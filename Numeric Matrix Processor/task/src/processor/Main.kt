package processor

import kotlin.math.pow
import kotlin.system.exitProcess

fun main() {
    NumericProcessor.run()
}

object NumericProcessor {
    fun run() {
        while (true) {
            println(
                """
                1. Add matrices
                2. Multiply matrix by a constant
                3. Multiply matrices
                4. Transpose matrix
                5. Calculate a determinant
                6. Inverse matrix
                0. Exit
            """.trimIndent()
            )

            print("Your choice: ")

            when (readln().toInt()) {
                1 -> addMatrices()
                2 -> multiplyByConstant()
                3 -> multiplyMatrices()
                4 -> transposeMatrix()
                5 -> calculateDeterminant()
                6 -> inverseMatrix()
                0 -> exitProcess(0)
            }
        }
    }

    private fun inverseMatrix() {
        print("Enter matrix size: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")
        
        if (matrix.determinant() == 0.0) {
            println("This matrix doesn't have an inverse.")
            return
        }
        
        println(matrix.inverse())
    }

    private fun calculateDeterminant() {
        print("Enter matrix size: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")

        val determinant = matrix.determinant()

        if (determinant % 1 > 0) {
            println(determinant)
        } else {
            println(determinant.toInt())
        }
    }

    private fun transposeMatrix() {
        println("""
            
            1. Main diagonal
            2. Side diagonal
            3. Vertical line
            4. Horizontal line
        """.trimIndent())

        println("Your choice: ")
        when (readln().toInt()) {
            1 -> transposeAlongMainDiagonal()
            2 -> transposeAlongSideDiagonal()
            3 -> transposeAlongVerticalLine()
            4 -> transposeAlongHorizontalLine()
        }
    }

    private fun transposeAlongSideDiagonal() {
        print("Enter size of matrix: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")
        println(matrix.transposeAlongSideDiagonal())
    }

    private fun transposeAlongMainDiagonal() {
        print("Enter size of matrix: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")
        println(matrix.transposeAlongMainDiagonal())
    }

    private fun transposeAlongHorizontalLine() {
        print("Enter size of matrix: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")
        println(matrix.transposeAlongHorizontalLine())
    }

    private fun transposeAlongVerticalLine() {
        print("Enter size of matrix: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix = createMatrix(rows, cols)

        println("The result is:")
        println(matrix.transposeAlongVerticalLine())
    }

    private fun multiplyMatrices() {
        print("Enter size of first matrix: ")
        val (rows1, cols1) = readln().split(" ").map { it.toInt() }
        println("Enter first matrix:")
        val matrix1 = createMatrix(rows1, cols1)

        print("Enter size of second matrix: ")
        val (rows2, cols2) = readln().split(" ").map { it.toInt() }
        println("Enter second matrix:")
        val matrix2 = createMatrix(rows2, cols2)

        if (!matrix1.canBeMultipliedBy(matrix2)) {
            println("The operation cannot be performed.\n")
            return
        }

        println("The result is:")
        println(matrix1 * matrix2)
        println()
    }

    private fun multiplyByConstant() {
        print("Enter size of matrix: ")
        val (rows, cols) = readln().split(" ").map { it.toInt() }
        println("Enter matrix:")
        val matrix1 = createMatrix(rows, cols)

        print("Enter constant: ")
        val scalar = readln().toDouble()

        println("The result is:")
        println(matrix1 * scalar)
        println()
    }

    private fun addMatrices() {
        print("Enter size of first matrix: ")
        val (rows1, cols1) = readln().split(" ").map { it.toInt() }
        println("Enter first matrix:")
        val matrix1 = createMatrix(rows1, cols1)

        print("Enter size of second matrix: ")
        val (rows2, cols2) = readln().split(" ").map { it.toInt() }
        println("Enter second matrix:")
        val matrix2 = createMatrix(rows2, cols2)

        if (matrix1.size() != matrix2.size()) {
            println("The operation cannot be performed.\n")
            return
        }

        println("The result is:")
        println(matrix1 + matrix2)
        println()
    }

    private fun createMatrix(rows: Int, cols: Int): Matrix {
        val matrix = Matrix(rows, cols)

        val content = mutableListOf<String>()

        repeat(rows) {
            content.add(readln())
        }

        matrix.populateMatrix(content)

        return matrix
    }
}


class Matrix(private var rows: Int, private var columns: Int) {
    private var matrix = MutableList(size = rows) { MutableList(size = columns) { 0.0 } }

    fun populateMatrix(content: MutableList<String>) {
        for ((i, row) in content.withIndex()) {
            for ((j, number) in row.split(" ").withIndex()) {
                matrix[i][j] = number.toDouble()
            }
        }
    }

    operator fun plus(other: Matrix): Matrix {
        val newMatrix = Matrix(rows, columns)

        for ((i, row) in matrix.withIndex()) {
            for ((j, column) in row.withIndex()) {
                newMatrix[i][j] = column + other[i][j]
            }
        }

        return newMatrix
    }

    operator fun get(i: Int): MutableList<Double> {
        return matrix[i]
    }

    override fun toString(): String {
        val result = StringBuilder()

        if (matrix.all { row -> row.all { it % 1 == 0.0 } }) {
            matrix.forEach { row -> result.appendLine(row.map { it.toInt() }.joinToString(" ")) }
        } else {
            matrix.forEach { result.appendLine(it.joinToString(" ")) }
        }

        return result.toString()
    }

    fun size(): Pair<Int, Int> {
        return Pair(rows, columns)
    }

    operator fun times(scalar: Int): Matrix {
        return this * scalar.toDouble()
    }

    operator fun times(scalar: Double): Matrix {
        val newMatrix = Matrix(rows, columns)

        for ((i, row) in matrix.withIndex()) {
            for ((j, column) in row.withIndex()) {
                newMatrix[i][j] = column * scalar
            }
        }

        return newMatrix
    }

    operator fun times(other: Matrix): Matrix {
        val newMatrix = Matrix(this.size().first, other.size().second)
        for (i in 0 until newMatrix.rows) {
            for (j in 0 until newMatrix.columns) {
                var element = 0.0
                for (x in 0 until this.columns) {
                    element += this[i][x] * other[x][j]
                }
                newMatrix[i][j] = element
            }
        }

        return newMatrix
    }

    fun canBeMultipliedBy(other: Matrix): Boolean {
        return this.size().second == other.size().first
    }

    fun transposeAlongMainDiagonal(): Matrix {
        val newMatrix = Matrix(rows = columns, columns = rows)

        for (i in 0 until rows) {
            for (j in 0 until columns) {
                newMatrix[j][i] = matrix[i][j]
            }
        }

        return newMatrix
    }

    fun transposeAlongSideDiagonal(): Matrix {
        return this.transposeAlongMainDiagonal().transposeAlongVerticalLine().transposeAlongHorizontalLine()
    }

    fun transposeAlongVerticalLine(): Matrix {
        val newMatrix = Matrix(rows, columns)

        for ((i, row) in matrix.withIndex()) {
            row.reverse()
            newMatrix[i] = row
        }

        return newMatrix
    }

    private operator fun set(i: Int, value: MutableList<Double>) {
        matrix[i] = value
    }

    fun transposeAlongHorizontalLine(): Matrix {
        val newMatrix = Matrix(rows, columns)

        newMatrix.matrix = matrix.toMutableList()

        newMatrix.matrix.reverse()

        return newMatrix
    }

    fun determinant(): Double {
        return Companion.determinant(this)
    }

    fun inverse(): Matrix {
        return this.adjoint() * (1 / determinant())
    }

    private fun adjoint(): Matrix {
        val adjoint = Matrix(rows, columns)

        for (i in 0 until rows) {
            for (j in 0 until columns) {
                adjoint[i][j] = getMinorAtIndex(copy(this), i, j).determinant() * (-1.0).pow(i + j)
            }
        }

        return adjoint.transposeAlongMainDiagonal()
    }

    companion object {
        private fun determinant(matrix: Matrix): Double {
            if (matrix.rows == 2 && matrix.columns == 2) {
                return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
            } else if (matrix.rows == 1 && matrix.columns == 1) {
                return matrix[0][0]
            }

            val minors = mutableListOf<Matrix>()

            for (i in 0 until matrix.columns) {
                minors.add(getMinorAtIndex(copy(matrix), 0, i))
            }

            return minors.mapIndexed { index, it -> matrix[0][index] * it.determinant() * (-1.0).pow(index) }.sum()
        }

        private fun getMinorAtIndex(matrix: Matrix, row: Int, column: Int): Matrix {
            matrix.matrix.removeAt(row)
            matrix.matrix.forEach { it.removeAt(column) }
            matrix.rows -= 1
            matrix.columns -= 1

            return matrix
        }

        private fun copy(matrix: Matrix): Matrix {
            val newMatrix = Matrix(matrix.rows, matrix.columns)

            for (row in 0 until matrix.rows) {
                for (col in 0 until matrix.columns) {
                    newMatrix[row][col] = matrix[row][col]
                }
            }

            return newMatrix
        }
    }
}