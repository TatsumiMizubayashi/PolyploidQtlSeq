namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Parent1 MostAlleleRateしきい値
    /// </summary>
    internal class Parent1MostAlleleRateThreshold
    {
        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最大値
        /// </summary>
        private const double MAXIMUM = 1.0;

        /// <summary>
        /// Parent1 MostAlleleRateしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Parent1MostAlleleRateThreshold(double threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// Parent1 MaxAlleleRateしきい値を取得する。
        /// </summary>
        internal double Value { get; }
    }
}
