namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 染色体長のしきい値
    /// </summary>
    internal class ChrSizeThreshold
    {
        /// <summary>
        /// しきい値の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// しきい値の最大値
        /// </summary>
        private const int MAXIMUM = 100_000_000;


        /// <summary>
        /// 染色体長のしきい値
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public ChrSizeThreshold(int threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// 染色体サイズのしきい値を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
