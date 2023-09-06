namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// 最大Bulk SNP-indexのしきい値
    /// </summary>
    internal class MaxBulkSnpIndexThreshold
    {
        /// <summary>
        /// 最大Bulk SNP-indexのしきい値の最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値の最大値
        /// </summary>
        private const double MAXIMUM = 1.0;


        /// <summary>
        /// 最大Bulk SNP-indexのしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public MaxBulkSnpIndexThreshold(double threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値を取得する。
        /// </summary>
        internal double Value { get; }

    }
}
