namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// リードの最低長
    /// </summary>
    internal class ReadLengthRequired
    {
        /// <summary>
        /// リード最低長の最小値
        /// </summary>
        private const int MINIMUM = 10;

        /// <summary>
        /// リード最低長の最大値
        /// </summary>
        private const int MAXIMUM = 150;

        /// <summary>
        /// リードの最低長を作成する。
        /// </summary>
        /// <param name="length">リード最低長</param>
        public ReadLengthRequired(int length)
        {
            if (length < MINIMUM || length > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(length));

            Value = length;
        }

        /// <summary>
        /// リードの最低長を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-l {Value}";
        }
    }
}
