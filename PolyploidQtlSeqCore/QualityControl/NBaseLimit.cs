namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// N塩基数の上限値
    /// </summary>
    internal class NBaseLimit
    {
        /// <summary>
        /// N塩基数の最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// N塩基数の最大値
        /// </summary>
        private const int MAXIMUM = 30;


        /// <summary>
        /// N塩基数の最大値を作成する。
        /// </summary>
        /// <param name="limit">N塩基最大数</param>
        public NBaseLimit(int limit)
        {
            if (limit < MINIMUM || limit > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(limit));

            Value = limit;
        }

        /// <summary>
        /// N塩基数の最大値を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-n {Value}";
        }
    }
}
