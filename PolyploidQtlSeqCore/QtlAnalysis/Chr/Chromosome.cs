namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 染色体情報
    /// </summary>
    internal class Chromosome
    {
        /// <summary>
        /// 染色体情報を作成する。
        /// </summary>
        /// <param name="name">名前</param>
        /// <param name="length">長さ</param>
        public Chromosome(string name, int length)
        {
            ArgumentException.ThrowIfNullOrEmpty(name);
            ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);

            Name = name;
            Length = length;
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// 染色体長を取得する。
        /// </summary>
        public int Length { get; }

        /*
         * 本来なら染色体名だけで良いが、染色体名の中に-が入っている場合、 
         * 染色体名のみ指定するとエラーになる。
         * そのため、染色体名:1-Endという表記で使用する。
         */

        /// <summary>
        /// bcftools用の位置情報に変換する。
        /// </summary>
        /// <returns>位置情報</returns>
        internal string ToRegion()
        {
            return $"{Name}:1-{Length}";
        }
    }
}
