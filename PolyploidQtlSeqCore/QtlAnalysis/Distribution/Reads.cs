namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 調査変異のリード配列
    /// </summary>
    internal class Reads
    {
        /// <summary>
        /// 調査変異のリード配列を作成する。
        /// </summary>
        /// <param name="gtValues">GT値配列</param>
        public Reads(int[] gtValues)
        {
            var altCount = gtValues.Sum();

            SnpIndex = (double)altCount / gtValues.Length;
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public double SnpIndex { get; }

    }
}
