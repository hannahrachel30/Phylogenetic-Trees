import java.io.*;
import java.util.*;

public class PhylogeneticTrees {
    // Function to print the OPT matrix with amino acid names
    public static void printMatrix(double[][] matrix, List<String> labels) {
        int m = matrix.length;
        int n = matrix[0].length;
        int w = labels.get(0).length() + 10;

        System.out.println();

        // Print column headers
        System.out.printf("%" + w + "s", " ");
        for (int i = 0; i < m; i++) {
            System.out.printf("%" + w + "s", labels.get(i));
        }
        System.out.println();

        // Print the OPT matrix with amino acid names
        for (int i = 0; i < m; i++) {
            System.out.printf("%" + w + "s", labels.get(i));
            for (int j = 0; j < n; j++) {
                System.out.printf("%" + w + ".2f", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    // Function to read the matrix from a file
    public static double[][] readMatrixFromFile(String filename, List<String> labels) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line;
        if ((line = reader.readLine()) != null) {
            int n = Integer.parseInt(line.trim());
            if ((line = reader.readLine()) != null) {
                String[] labelArr = line.split("\t");
                labels.addAll(Arrays.asList(labelArr)); // Extracting labels from the first line of the file
            }
            double[][] matrix = new double[n][n]; // Resize the matrix to n x n
            for (int i = 0; i < n; i++) {
                if ((line = reader.readLine()) != null) {
                    String[] values = line.split("\t");
                    for (int j = 0; j < n; j++) {
                        matrix[i][j] = Double.parseDouble(values[j]); // Populating the matrix with values from the file
                    }
                }
            }
            reader.close();
            return matrix;
        }
        reader.close();
        return new double[0][0];
    }

    // Function to find the indices of the minimum value in the distance matrix
    public static void findMinimum(double[][] transition, int[] minIndices, double[] minVal, double[][] matrix, double[] rValues) {
        minVal[0] = 9999.99; // Setting an initial large value for minimum comparison
        minIndices[0] = 0;
        minIndices[1] = 0;
        int n = transition.length;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double temp = (matrix[i][j] + rValues[i] - rValues[j]);
                double minTemp = (matrix[minIndices[0]][minIndices[1]] + rValues[minIndices[0]] - rValues[minIndices[1]]);
                if (transition[i][j] < minVal[0]) {
                    minVal[0] = transition[i][j]; // Update minimum value and corresponding indices
                    minIndices[0] = i;
                    minIndices[1] = j;
                } else if (transition[i][j] == minVal[0] && temp < minTemp) {
                    minVal[0] = transition[i][j]; // Update minimum value and corresponding indices
                    minIndices[0] = i;
                    minIndices[1] = j;
                }
            }
        }
    }

    // Function implementing the Neighbor Joining algorithm
    public static void neighborJoin(double[][] matrix, List<String> labels) {
        List<Integer> clusterItems = new ArrayList<>(Collections.nCopies(labels.size(), 1)); // Track the number of items in each cluster

        int n = matrix.length;

        while (n > 2) {
            double[] rValues = new double[n];

            System.out.println("r-values:");
            // Compute r-values
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    rValues[i] += matrix[i][j];
                }
                rValues[i] /= (n - 2);
                System.out.printf("%s : %.2f\t", labels.get(i), rValues[i]);
            }
            System.out.println();

            // compute transition matrix
            double[][] transition = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    transition[i][j] = matrix[i][j] - rValues[i] - rValues[j];
                    transition[j][i] = transition[i][j];
                }
            }

            System.out.println("Transition matrix:");
            printMatrix(transition, labels);

            n--;

            double[] minimum = {9999};
            int[] minIndices = {0, 0};

            findMinimum(transition, minIndices, minimum, matrix, rValues); // Find the minimum distance pair

            int i = minIndices[0], j = minIndices[1];
            System.out.printf("Clusters to merge - %s and %s\n", labels.get(i), labels.get(j));
            System.out.printf("Distances to ancestor - %s : %.2f, %s : %.2f\n",
                    labels.get(i), (matrix[i][j] + rValues[i] - rValues[j]) / 2,
                    labels.get(j), (matrix[i][j] + rValues[j] - rValues[i]) / 2);

            String newCluster = "(" + labels.get(i) + ", " + labels.get(j) + ")"; // Create a new cluster name
            List<String> newLabels = new ArrayList<>();
            double[][] newMatrix = new double[n][n]; // Resize the new matrix

            newLabels.add(newCluster); // Add the new cluster to the labels

            int n1 = clusterItems.get(i);
            int n2 = clusterItems.get(j);

            // Update the number of items in the clusters
            clusterItems.remove(Math.max(i, j));
            clusterItems.remove(Math.min(i, j));
            clusterItems.add(0, n1 + n2);

            // Reconstruct the new label list and new distance matrix
            for (int k = 0; k < labels.size(); k++) {
                if (k != i && k != j) {
                    newLabels.add(labels.get(k));
                }
            }

            int c1 = 1;
            // Update the first row and column of the new matrix
            for (int m = 0; m < matrix.length; m++) {
                if (m != i && m != j) {
                    newMatrix[c1][0] = (matrix[m][i] + matrix[m][j] - matrix[i][j]) / 2;
                    newMatrix[0][c1] = newMatrix[c1][0];
                    c1++;
                }
            }

            // Update the remaining cells of the new matrix
            int c2 = 2;
            c1 = 1;
            for (int m = 0; m < matrix.length; m++) {
                if (m != i && m != j) {
                    for (int p = m + 1; p < matrix.length; p++) {
                        if (p != i && p != j) {
                            newMatrix[c2][c1] = matrix[m][p];
                            newMatrix[c1][c2] = newMatrix[c2][c1];
                            c2++;
                        }
                    }
                    c1++;
                    c2 = c1 + 1;
                }
            }

            matrix = newMatrix; // Update the distance matrix
            labels = newLabels; // Update the labels with the new cluster

            System.out.println("Distance Matrix:");
            printMatrix(matrix, labels); // Print the current distance matrix
        }

        System.out.println("Distance between the remaining clusters:");
        System.out.printf("%s : %.2f, %s : %.2f\n", labels.get(0), matrix[0][1], labels.get(1), matrix[0][1]);

        // Print the final tree
        System.out.printf("Resulting Newick format - %s %s\n", labels.get(0), labels.get(1));
    }

    public static void main(String[] args) throws IOException {
        List<String> labels = new ArrayList<>();

        System.out.print("Enter the name of the distance matrix file: ");
        Scanner scanner = new Scanner(System.in);
        String filename = scanner.nextLine(); // Input the filename

        double[][] matrix = readMatrixFromFile(filename, labels); // Read the matrix from the file

        if (matrix.length == 0) {
            System.out.println("Failed to read the matrix from the file.");
            return;
        }

        System.out.print("Enter the method to run (UPGMA or NJ): ");
        String method = scanner.nextLine();

        if (method.equals("UPGMA")) {
            // Running UPGMA on the provided distance matrix
            System.out.println("Running UPGMA on the matrix:");
            printMatrix(matrix, labels); // Print the initial distance matrix
            UPGMA(matrix, labels); // Perform the UPGMA algorithm
        } else if (method.equals("NJ")) {
            // Running Neighbour Join on the provided distance matrix
            System.out.println("Running Neighbor Join on the matrix:");
            printMatrix(matrix, labels); // Print the initial distance matrix
            neighborJoin(matrix, labels); // Perform the Neighbor Join algorithm
        }

        scanner.close();
    }

    // Function implementing the UPGMA algorithm (unchanged)
    public static void UPGMA(double[][] matrix, List<String> labels) {
        List<Integer> clusterItems = new ArrayList<>(Collections.nCopies(labels.size(), 1)); // Track the number of items in each cluster

        while (labels.size() > 1) { // Continue until only one cluster remains
            double minimum = 9999;
            int i = 0, j = 0;
            int n = matrix.length;

            int[] minIndices = {i, j};
            double[] minVal = {minimum};
            findMinimum(matrix, minIndices, minVal, matrix, new double[matrix.length]); // Find the minimum distance pair
            i = minIndices[0];
            j = minIndices[1];

            // Print the clusters being merged and their distance
            System.out.printf("Clusters to merge - %s and %s with distance = %.2f\n", labels.get(i), labels.get(j), matrix[i][j]);

            String newCluster = "(" + labels.get(i) + ", " + labels.get(j) + ")"; // Create a new cluster name
            List<String> newLabels = new ArrayList<>();
            double[][] newMatrix = new double[matrix.length - 1][matrix.length - 1]; // Resize the new matrix

            newLabels.add(newCluster); // Add the new cluster to the labels

            int n1 = clusterItems.get(i);
            int n2 = clusterItems.get(j);

            // Update the number of items in the clusters
            if (i > j) {
                clusterItems.remove(i);
                clusterItems.remove(j);
            } else {
                clusterItems.remove(j);
                clusterItems.remove(i);
            }

            clusterItems.add(0, n1 + n2);

            // Reconstruct the new label list and new distance matrix
            for (int k = 0; k < labels.size(); k++) {
                if (k != i && k != j) {
                    newLabels.add(labels.get(k));
                }
            }

            int c1 = 1;
            // Update the first row and column of the new matrix
            for (int m = 0; m < matrix.length; m++) {
                if (m != i && m != j) {
                    newMatrix[c1][0] = (matrix[m][i] * n1 + matrix[m][j] * n2) / (n1 + n2);
                    newMatrix[0][c1] = newMatrix[c1][0];
                    c1++;
                }
            }

            // Update the remaining cells of the new matrix
            int c2 = 2;
            c1 = 1;
            for (int m = 0; m < matrix.length; m++) {
                if (m != i && m != j) {
                    for (int p = m + 1; p < matrix.length; p++) {
                        if (p != i && p != j) {
                            newMatrix[c2][c1] = matrix[m][p];
                            newMatrix[c1][c2] = newMatrix[c2][c1];
                            c2++;
                        }
                    }
                    c1++;
                    c2 = c1 + 1;
                }
            }

            matrix = newMatrix; // Update the distance matrix
            labels = newLabels; // Update the labels with the new cluster

            printMatrix(matrix, labels); // Print the current distance matrix
        }

        // Print the final clustering result in Newick format
        System.out.printf("Resulting Newick format - %s\n", labels.get(0));
    }
}
