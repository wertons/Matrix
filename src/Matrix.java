public class Matrix {

    static double trace(double[][] mat) {
        try {
            checkMatrix(mat);
            testSquare(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Find the lowest dimension of the matrix, then add all elements of the main diagonal
        int lowestLength = 0;
        double result = 0.0;
        lowestLength = Math.min(mat.length, mat[0].length);
        for (int i = 0; i < lowestLength; i++) {
            result += mat[i][i];
        }

        return result;
    }

    static double[][] add(double[][] mat1, double[][] mat2) {
        try {
            checkMatrix(mat1);
            checkMatrix(mat2);
            equalMatrix(mat1, mat2);

        } catch (Exception e) {
            e.printStackTrace();
        }
        //Create a new array copy of the first one to get the required size, as the addition needs to be performed
        // on two equal arrays, then simply do basic additions on a 2 deep loop
        for (int i = 0; i < mat1.length; i++) {
            for (int j = 0; j < mat1[i].length; j++) {
                mat1[i][j] = mat1[i][j] + mat2[i][j];
            }
        }
        return mat1;
    }

    static double[][] mult(double[][] mat1, double[][] mat2) {
        try {
            checkMatrix(mat1);
            checkMatrix(mat2);
            testMult(mat1, mat2);

        } catch (Exception e) {
            e.printStackTrace();
        }
        //We first create the result array that is the width of mat1 and length of mat2
        double[][] matresult = new double[mat1.length][mat2[0].length];

        for (int i = 0; i < matresult.length; i++) {
            for (int j = 0; j < matresult[0].length; j++) {
                //We call a external method
                matresult[i][j] = multiplyMatrixCells(mat1, mat2, i, j);
            }
        }
        return matresult;
    }

    private static double multiplyMatrixCells(double[][] mat1, double[][] mat2, int row, int col) {
        try {
            checkMatrix(mat1);
            checkMatrix(mat2);
            equalMatrix(mat1, mat2);

        } catch (Exception e) {
            e.printStackTrace();
        }
        double cell = 0.0;
        //We find the minimum amount of operations for the multiplication by finding the shortest between the width
        // of mat1 and length of mat2
        double lengthiest = 0;
        lengthiest = Math.max(mat1[0].length, mat2.length);
        //We add the multiplications of the rows and columns to the result
        for (int i = 0; i < lengthiest; i++) {
            cell += mat1[row][i] * mat2[i][col];
        }
        return cell;
    }

    static double[][] mult(double[][] mat, double n) {
        try {
            checkMatrix(mat);

        } catch (Exception e) {
            e.printStackTrace();
        }
        //We multiply all cells by the multiplier
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                mat[i][j] = mat[i][j] * n;
            }
        }
        return mat;
    }

    static double[][] power(double[][] mat, int p) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //We simply copy the array and multiply on the copy n times
        double[][] result = mat;

        if (p == 0) {
            result = identityMatrix(mat);
        } else {
            for (int i = 0; i < p - 1; i++) {
                result = mult(result, mat);
            }
        }
        return result;
    }

    static double[][] div(double[][] mat1, double[][] mat2) {
        try {
            checkMatrix(mat1);
            checkMatrix(mat2);
            equalMatrix(mat1, invert(mat2));
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (invert(mat2) == null) {
            return null;
        }
        return mult(mat1, invert(mat2));
    }

    static double[][] submatrix(double[][] mat, int x1, int y1, int x2, int y2) {
        try {
            checkMatrix(mat);

        } catch (Exception e) {
            e.printStackTrace();
        }
        if (x1 < 0 || y1 < 0 || x2 > mat.length - 1 || y2 > mat[0].length - 1) {
            try {
                throw new Exception("Invalid sub-matrix: The positions of the sub-matrix have to be in the range of the original matrix");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        //First create the array based on the new size
        double[][] result = new double[y2 - y1 + 1][x2 - x1 + 1];

        int resX = 0;
        int resY = 0;
        boolean incre = false;
        for (int i = 0; i < mat[0].length; i++) {
            for (int j = 0; j < mat.length; j++) {
                //Check if the position is the required one to extract, but make a independent counter for the second array
                if (i >= x1 && i <= x2 && j >= y1 && j <= y2) {
                    result[resX][resY] = mat[j][i];
                    resX++;
                    incre = true;
                }
            }
            //For resets and increments of the seconds set of incremental base it only when
            resX = 0;
            if (incre) {
                resY++;
            }
            incre = false;
        }
        return result;
    }

    static double[][] invert(double[][] mat) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][] minors = calcMinors(mat);
        double[][] result = new double[mat.length][mat.length];
        //Apply the logic of +'s and -'s tp the cells
        for (int i = 0; i < mat.length; i++) {
            System.arraycopy(minors[i], 0, result[i], 0, mat.length);

        }
        //We first transpose the result
        result = transpose(result);
        //We obtain the dividend via the 1 / the determinant of the original matrix
        try {
            testDeterminant(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double div = 1.0 / determinant(mat);
        //We multiply the dividend by the transposed matrix
        mult(result, div);
        //Apply the logic of +'s and -'s tp the cells
        int counter = 0;
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                result[i][j] = (Math.pow(-1, counter)) * result[i][j];
                counter++;
            }
            if (i % 2 == 0) {
                counter = 1;
            } else {
                counter = 0;
            }
        }
        if (div == Double.POSITIVE_INFINITY) {
            return null;
        } else {
            return result;
        }
    }

    static double[][] identityMatrix(double[][] mat) {
        try {
            checkMatrix(mat);
            testSquare(mat);

        } catch (Exception e) {
            e.printStackTrace();
        }
        //Construct a identity matrix, basically a matrix of same size but with the main diagonal being all 1's and
        // the rest 0`s
        double[][] result = new double[mat.length][mat.length];
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                if (i == j) {
                    result[i][j] = 1;
                } else {
                    result[i][j] = 0;
                }
            }
        }
        return result;
    }

    static double[][] calcMinors(double[][] mat) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //Create a matrix of minors, which is basically calling the getMinor function once for each cell, then
        // calculating the determinant of each minor
        double[][] result = new double[mat.length][mat[0].length];
        double a = 0.0;
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                try {
                    testDeterminant(getMinor(mat, i, j));
                } catch (Exception e) {
                    e.printStackTrace();
                }
                result[i][j] = determinant(getMinor(mat, i, j));
            }
        }
        return result;
    }

    static double[][] getMinor(double[][] mat, int x, int y) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (x < 0 || y < 0 || x > mat.length - 1 || y > mat[0].length - 1) {
            try {
                throw new Exception("Invalid sub-matrix: The positions of the sub-matrix have to be in the range of the original matrix");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        //We create a new array with size -1 as related to the original
        double[][] result = new double[mat.length - 1][mat[0].length - 1];

        //Same logic as a submatrix
        int resX = 0;
        int resY = 0;
        boolean incre = false;

        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                if (i != x && j != y) {
                    result[resX][resY] = mat[i][j];
                    resY++;
                    incre = true;
                }
            }
            //For resets and increments of the seconds set of incremental base it only when
            resY = 0;
            if (incre) {
                resX++;
            }
            incre = false;
        }
        return result;
    }

    static double determinant(double[][] mat) {
        try {
            checkMatrix(mat);
            testSquare(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double result = 0.0;
        //If length 1 can't be calculated
        if (mat.length == 1) {
            result = mat[0][0];
            //If length 2 just call the determinant
        } else if (mat.length == 2) {
            result = subDeterminant(mat);
        } else {
            //Calc determinant by calculating the minors then solving
            for (int i = 0; i < mat.length; i++) {
                double[][] tmpMat = getMinor(mat, 0, i);
                try {
                    testDeterminant(tmpMat);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                result += mat[0][i] * determinant(tmpMat) * Math.pow(-1, i);
            }
        }
        return result;
    }

    static double subDeterminant(double[][] mat) {
        try {
            checkMatrix(mat);
            testSquare(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //We calculate a 2 by 2 determinant
        return (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]);
    }

    static double[][] transpose(double[][] mat) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //We create the transposed array by flipping the lengths, in case of a non-equilateral array
        double[][] matResult = new double[mat[0].length][mat.length];

        //We flip the arrays by switching the array replacement order: [i][j] -> [j][i]
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                matResult[i][j] = mat[j][i];
            }
        }
        return matResult;
    }

    static boolean isOrtho(double[][] mat) {
        try {
            checkMatrix(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //A matrix is orthogonal if the transposed matrix and the inverse matrix are equal, so we check if the
        // transpose and the inversion are equal, we also check if the determinant is 0 in which case this operation
        // is impossible
        boolean result = true;
        try {
            testDeterminant(mat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][] transMat = transpose(mat);
        double[][] invMat = invert(mat);
        negativeZeros(transMat);
        negativeZeros(invMat);
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                assert invMat != null;
                if (transMat[i][j] != invMat[i][j] ) {
                    result = false;
                    break;
                }
            }
        }
        return result;
    }

    static double[] cramer(double[][] mat) {
        try {
            checkMatrix(mat);
            if (mat.length != mat[0].length -1){
                throw new Exception("Non Cramer logic matrix: The matrix has to be a square with and additional column in order to apply Cramer's Rule");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        double[][] subMat = submatrix(mat, 0, 0, mat.length - 1, mat[0].length - 2);
        try {
            testDeterminant(subMat);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double mainDet = determinant(subMat);
        double[] result = new double[mat.length];
        double[][] tmpSubMat = subMat;

        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                tmpSubMat[j][i] = mat[j][mat[0].length - 1];
            }
            result[i] = determinant(tmpSubMat) / mainDet;
            tmpSubMat = submatrix(mat, 0, 0, mat.length - 1, mat[0].length - 2);
        }

        return result;
    }

    // Funció que mostra una matriu per pantalla
    static void printMat(double[][] mat) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                System.out.printf("%06.2f ", mat[i][j]);
            }
            System.out.println();
        }
    }

    static void checkMatrix(double[][] mat) throws Exception {
        //Checks if there is an actual matrix bigger than 0x0
        if (mat.length <= 0) {
            throw new Exception("Null size: Input has to be a matrix bigger than 0");
        }

        //Checks if the input is an array by comparing all columns
        for (int i = 0; i < mat.length; i++) {
            if (mat.length <= 0 || mat[0].length <= 0) {
                throw new Exception("Null size: Input has to be a matrix bigger than 0");
            }

            if (mat[i].length != mat[0].length) {
                throw new Exception("Non-matrix input: The input has to be a rectangular matrix");
            }
        }
    }

    static void testSquare(double[][] mat) throws Exception {
        for (int i = 0; i < mat.length; i++) {
            if (mat[0].length != mat.length) {
                throw new Exception("Non-square matrix: The matrix is not a square");
            }
        }
    }

    static void equalMatrix(double[][] mat1, double[][] mat2) throws Exception {
        if (mat1.length != mat2.length)
            for (int i = 0; i < mat1.length; i++) {
                if (mat1[i].length != mat2[i].length) {
                    throw new Exception("Non-equal matrices: The matrices are not equal");
                }
            }
    }
    static void testMult(double[][] mat1, double[][] mat2) throws Exception {
                if (mat1.length != mat2[0].length) {
                    throw new Exception("Non-operable matrices: The matrices are right for realizing a multiplication");
                }

    }
    static double[][] negativeZeros(double[][] mat){
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                if (mat[i][j] == -0.0){
                    mat[i][j] = 0;
                }
            }
        }
        return mat;
    }

    static void testDeterminant(double[][] mat) throws Exception {
        //Calculate the determinant
        double det = determinant(mat);

        //Checks if the determinant is either null or infite, in which cases operations are imposible
        if (det == Double.POSITIVE_INFINITY || Double.isNaN(det)) {
            throw new Exception("Inconsistent determinant: The determinant is either non-existent or has a infinite value");
        }

    }

}