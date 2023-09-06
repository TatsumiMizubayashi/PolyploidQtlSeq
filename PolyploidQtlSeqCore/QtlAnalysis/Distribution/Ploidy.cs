namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 倍数性
    /// </summary>
    internal class Ploidy
    {
        /// <summary>
        /// ploidyの最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// ploidyの最大値
        /// </summary>
        private const int MAXIMUM = 20;


        /// <summary>
        /// 倍数性を作成する。
        /// </summary>
        /// <param name="ploidy">倍数性</param>
        public Ploidy(int ploidy)
        {
            if (ploidy < MINIMUM || ploidy > MAXIMUM) throw new ArgumentException(null, nameof(ploidy));
            if (ploidy % 2 != 0) throw new ArgumentException(null, nameof(ploidy));

            Value = ploidy;
        }

        /// <summary>
        /// 倍数性を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
