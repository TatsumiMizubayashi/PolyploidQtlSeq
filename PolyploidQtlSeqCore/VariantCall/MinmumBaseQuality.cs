namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// BaseQuality最低値
    /// </summary>
    internal class MinmumBaseQuality
    {
        /// <summary>
        /// min-BQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// min-BQの最大値
        /// </summary>
        private const int MAXIMUM = 60;

        /// <summary>
        /// Base Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minBq">Base Quality最小値</param>
        public MinmumBaseQuality(int minBq)
        {
            if (minBq < MINIMUM || minBq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(minBq));

            Value = minBq;
        }

        /// <summary>
        /// Base Quality最小値を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
